from __future__ import absolute_import, unicode_literals

import logging
import pickle
from datetime import timedelta
from typing import List, Optional

import pandas as pd
from billiard.exceptions import SoftTimeLimitExceeded
from celery import shared_task, Task
from cleanup_later.models import CleanupFile
from django.core.mail import EmailMultiAlternatives
from django.db.models import F

from viewer.models import User, EnrichmentResult
from viewer.src.constants import make_gsea_export_directories, GENE_SYMBOL
from viewer.src.gsea import perform_gsea, perform_prerank
from viewer.src.ora import run_ora
from viewer.src.utils import read_data_file


@shared_task
def email(subject: str, text_content: str, html_content: str, sender: str, recipient_list: List[str]):
    msg = EmailMultiAlternatives(subject, text_content, sender, recipient_list)
    msg.attach_alternative(html_content, "text/html")
    msg.send()
    return None


@shared_task(soft_time_limit=28800)
def deploy_gsea(
    data_path: str,
    class_labels_path: str,
    data_filename: str,
    class_filename: str,
    gmt_files: list,
    output_dir: str,
    min_size: int,
    max_size: int,
    method: str,
    permutation_type: str,
    permutation_num: int,
    user_mail: str,
    job_id: str,
    read_counts_path: Optional = None,
    read_counts_filename: Optional = None,
):
    current_user = User.objects.filter(email=user_mail)[0]
    job = EnrichmentResult.objects.filter(result_id=job_id)[0]
    err = None

    try:
        make_gsea_export_directories()

        # Check if quota has been reached
        if not current_user.is_quota_left():
            job.result_status = 0
            job.error_message = "You have exceeded your quota for the number of experiments that can be run. " \
                                "Please cancel an existing experiment to start a new one."

            job.save()

            return False

        # Remove transient files
        CleanupFile.register(data_path, timedelta(minutes=30))
        CleanupFile.register(class_labels_path, timedelta(minutes=30))

        logging.info(f'Loading data file {data_filename}....')

        try:
            df = read_data_file(data_path, data_filename)
        except:
            raise ValueError(open(data_path, "r").read())

        if GENE_SYMBOL in df:
            df.set_index(GENE_SYMBOL, inplace=True)

        logging.info(f'Loading class labels file {class_filename}....')

        try:
            class_df = read_data_file(class_labels_path, class_filename)
        except:
            raise ValueError(open(class_labels_path, "r").read())

        logging.info('Files to run GSEA have been loaded....')

        # Get class labels vector
        class_vector = class_df['class_label'].to_list()

        # Run DESeq2
        if read_counts_path:

            CleanupFile.register(read_counts_path, timedelta(minutes=30))

            logging.info(f'Loading read counts file {read_counts_filename}....')

            try:
                read_counts_df = read_data_file(read_counts_path, read_counts_filename)
            except:
                raise ValueError(open(read_counts_path, "r").read())

            logging.info('Read counts file has been loaded....')

            if read_counts_df.index.name == GENE_SYMBOL:
                read_counts_df.reset_index(inplace=True)

            deseq = DeseqTask(
                count_matrix=read_counts_df,
                design_matrix=class_df,
                design_formula='~ class_label',
                gene_column=GENE_SYMBOL
            )

            pd_from_r_df = deseq.run()

            job.fold_change_results = pickle.dumps(pd_from_r_df)

        results = []

        for gmt_file in gmt_files:
            logging.info(f'Running GSEA on gene sets from {gmt_file}')

            gsea_results = perform_gsea(
                data=df,
                gmt=gmt_file,
                class_vector=class_vector,
                output_dir=output_dir,
                min_size=min_size,
                max_size=max_size,
                method=method,
                permutation_type=permutation_type,
                permutation_num=permutation_num,
            )

            logging.info('Getting results...')

            results_df = gsea_results.res2d
            results.append(results_df)

        df = pd.concat(results)

        job.result = pickle.dumps(df)
        job.result_status = 2

    except SoftTimeLimitExceeded:
        job.result_status = 0
        job.error_message = "The experiment exceeded the acceptable time limit (8 hours)"

    except Exception as e:
        err = e
        job.result_status = 0
        job.error_message = str(e)

    finally:
        User.objects.filter(email__exact=user_mail).update(num_of_jobs=F('num_of_jobs') - 1)
        job.save()

    return err


@shared_task(soft_time_limit=28800)
def deploy_prerank(
    rnk_path: str,
    rnk_filename: str,
    gmt_files: list,
    output_dir: str,
    min_size: int,
    max_size: int,
    permutation_num: int,
    user_mail: str,
    job_id: str,
    read_counts_path: Optional = None,
    read_counts_filename: Optional = None,
    class_labels_path: Optional = None,
    class_filename: Optional = None,
):
    current_user = User.objects.filter(email=user_mail)[0]
    job = EnrichmentResult.objects.filter(result_id=job_id)[0]
    err = None

    try:
        make_gsea_export_directories()

        # Check if quota has been reached
        if not current_user.is_quota_left():
            job.result_status = 0
            job.error_message = "You have exceeded your quota for the number of experiments that can be run. " \
                                "Please cancel an existing experiment to start a new one."

            job.save()

            return False

        # Remove transient files
        CleanupFile.register(rnk_path, timedelta(minutes=30))

        logging.info(f'Loading preranked file {rnk_path}....')

        try:
            df = read_data_file(rnk_path, rnk_filename)
        except:
            raise ValueError(open(rnk_path, "r").read())

        logging.info('Pre-ranked file has been loaded....')

        # Run DESeq2
        if read_counts_path:

            CleanupFile.register(read_counts_path, timedelta(minutes=30))
            CleanupFile.register(class_labels_path, timedelta(minutes=30))

            logging.info(f'Loading read counts file {read_counts_filename}....')

            try:
                read_counts_df = read_data_file(read_counts_path, read_counts_filename)
            except:
                raise ValueError(open(read_counts_path, "r").read())

            logging.info('Read counts file has been loaded....')

            logging.info(f'Loading class labels file {class_filename}....')

            try:
                class_df = read_data_file(class_labels_path, class_filename)
            except:
                raise ValueError(open(class_labels_path, "r").read())

            logging.info('Class labels file has been loaded....')

            if read_counts_df.index.name == GENE_SYMBOL:
                read_counts_df.reset_index(inplace=True)

            deseq = DeseqTask(
                count_matrix=read_counts_df,
                design_matrix=class_df,
                design_formula='~ class_label',
                gene_column=GENE_SYMBOL
            )

            pd_from_r_df = deseq.run()

            job.fold_change_results = pickle.dumps(pd_from_r_df)

        results = []

        for gmt_file in gmt_files:
            logging.info(f'Running GSEA Pre-Ranked on gene sets from {gmt_file}')

            prerank_results = perform_prerank(
                rnk=df,
                gmt=gmt_file,
                output_dir=output_dir,
                min_size=min_size,
                max_size=max_size,
                permutation_num=permutation_num,
            )

            logging.info('Getting results...')

            results_df = prerank_results.res2d

            results.append(results_df)

        df = pd.concat(results)

        job.result = pickle.dumps(df)
        job.result_status = 2

    except SoftTimeLimitExceeded:
        job.result_status = 0
        job.error_message = "The experiment exceeded the acceptable time limit (8 hours)"

    except Exception as e:
        err = e
        job.result_status = 0
        job.error_message = str(e)

    finally:
        User.objects.filter(email__exact=user_mail).update(num_of_jobs=F('num_of_jobs') - 1)
        job.save()

    return err


@shared_task(soft_time_limit=28800)
def deploy_ora(
    gmt_file_path: str,
    min_size: int,
    max_size: int,
    user_mail: str,
    job_id: str,
    set_gene_symbols: Optional[List[str]] = None,
    sig_threshold_fc: Optional = None,
    read_counts_path: Optional = None,
    design_matrix_path: Optional = None,
    fold_changes_path: Optional = None,
    read_counts_filename: Optional = None,
    design_matrix_filename: Optional = None,
    fold_changes_filename: Optional = None,

):
    current_user = User.objects.filter(email=user_mail)[0]
    job = EnrichmentResult.objects.filter(result_id=job_id)[0]
    err = None

    try:
        if not current_user.is_quota_left():
            job.result_status = 0
            job.error_message = "You have exceeded your quota for the number of experiments. Please cancel an " \
                                "existing experiment to start a new one."
            job.save()

            return False

        # Perform DGE analysis and run ORA on genes that pass significance
        if read_counts_path:

            # Remove transient files
            CleanupFile.register(read_counts_path, timedelta(minutes=30))
            CleanupFile.register(design_matrix_path, timedelta(minutes=30))

            logging.info(f'Loading data file {read_counts_filename}....')

            try:
                read_counts_df = read_data_file(read_counts_path, read_counts_filename)
            except:
                raise ValueError(open(read_counts_path, "r").read())

            if read_counts_df.index.name == GENE_SYMBOL:
                read_counts_df.reset_index(inplace=True)

            logging.info(f'Loading class labels file {design_matrix_filename}....')

            try:
                design_matrix = read_data_file(design_matrix_path, design_matrix_path)
            except:
                raise ValueError(open(read_counts_path, "r").read())

            logging.info('Files have been loaded. Starting DESeq2...')

            # Run DESeq2
            deseq = DeseqTask(
                read_counts_df,
                design_matrix,
                '~ class_label',
                GENE_SYMBOL
            )

            pd_from_r_df = deseq.run()

            logging.info("Filtering DEGs by adjusted p-value...")

            # Query genes which pass significance threshold based on adjusted p-value
            query_exp = f'{"padj"} <= {sig_threshold_fc}'

            df = pd_from_r_df.query(query_exp)

            ora_gene_set = set(df[GENE_SYMBOL])

            job.fold_change_results = pickle.dumps(pd_from_r_df)

            logging.info("Running ORA...")

            # Run ORA
            ora_results_df = run_ora(
                gmt_path=gmt_file_path,
                set_gene_symbols=ora_gene_set,
                min_size=min_size,
                max_size=max_size,
            )

            logging.info("ORA successfully run...")

            job.result = pickle.dumps(ora_results_df)

            # Job Success
            job.result_status = 2

        # Run ORA on DEGs that pass significance
        elif fold_changes_path:

            # Remove transient files
            CleanupFile.register(fold_changes_path, timedelta(minutes=30))

            logging.info(f'Loading fold changes file {fold_changes_filename}....')

            try:
                fold_changes_df = read_data_file(fold_changes_path, fold_changes_filename)
            except:
                raise ValueError(open(fold_changes_path, "r").read())

            logging.info(f'Fold changes file has been loaded....')

            logging.info("Filtering DEGs by adjusted p-value")

            # Query genes which pass significance threshold based on adjusted p-value
            query_exp = f'{"q_value"} <= {sig_threshold_fc}'

            df = fold_changes_df.query(query_exp)
            ora_gene_set = set(df[GENE_SYMBOL])

            logging.info("Running ORA")

            # Run ORA
            ora_results_df = run_ora(
                gmt_path=gmt_file_path,
                set_gene_symbols=ora_gene_set,
                min_size=min_size,
                max_size=max_size,
            )

            logging.info("ORA successfully run")

            job.result = pickle.dumps(ora_results_df)

            # Job Success
            job.result_status = 2

        # Run ORA
        else:
            ora_results_df = run_ora(
                gmt_path=gmt_file_path,
                set_gene_symbols=set(set_gene_symbols),
                min_size=min_size,
                max_size=max_size,
            )

            job.result = pickle.dumps(ora_results_df)

            # Job Success
            job.result_status = 2

    except SoftTimeLimitExceeded:
        job.result_status = 0
        job.error_message = "The experiment exceeded the acceptable time limit (8 hours)"

    except Exception as e:
        err = e
        job.result_status = 0
        job.error_message = str(e)

    finally:
        User.objects.filter(email__exact=user_mail).update(num_of_jobs=F('num_of_jobs') - 1)
        job.save()

    return err


@shared_task(soft_time_limit=28800)
def deploy_deseq(
    read_counts_path: str,
    read_counts_filename: str,
    design_matrix_path: str,
    design_matrix_filename: str,
    user_mail: str,
    job_id: str,
):
    current_user = User.objects.filter(email=user_mail)[0]
    job = EnrichmentResult.objects.filter(result_id=job_id)[0]
    err = None

    try:
        if not current_user.is_quota_left():
            job.result_status = 0
            job.error_message = "You have exceeded your quota for the number of experiments. Please cancel an " \
                                "existing experiment to start a new one."
            job.save()

            return False

        # Remove transient files
        CleanupFile.register(read_counts_path, timedelta(minutes=30))
        CleanupFile.register(design_matrix_path, timedelta(minutes=30))

        logging.info(f'Loading read counts file {read_counts_filename}....')

        try:
            read_counts_df = read_data_file(read_counts_path, read_counts_filename)
        except:
            raise ValueError(open(read_counts_path, "r").read())

        if read_counts_df.index.name == GENE_SYMBOL:
            read_counts_df.reset_index(inplace=True)

        logging.info(f'Loading design matrix file {design_matrix_path}....')

        try:
            design_matrix = read_data_file(design_matrix_path, design_matrix_filename)
        except:
            raise ValueError(open(design_matrix_path, "r").read())

        logging.info('Files have been loaded. Starting DESeq2...')

        deseq = DeseqTask(
            count_matrix=read_counts_df,
            design_matrix=design_matrix,
            design_formula='~ class_label',
            gene_column=GENE_SYMBOL
        )

        pd_from_r_df = deseq.run()

        job.fold_change_results = pickle.dumps(pd_from_r_df)
        job.result_status = 2

    except SoftTimeLimitExceeded:
        job.result_status = 0
        job.error_message = "The experiment exceeded the acceptable time limit (8 hours)"

    except Exception as e:
        err = e
        job.result_status = 0
        job.error_message = str(e)

    finally:
        User.objects.filter(email__exact=user_mail).update(num_of_jobs=F('num_of_jobs') - 1)
        job.save()

    return err


class DeseqTask(Task):
    def __init__(self, count_matrix, design_matrix, design_formula, gene_column):
        try:
            assert gene_column in count_matrix.columns, 'Wrong gene id column name'

        except AttributeError:
            raise Exception('Wrong Pandas dataframe?')

        import rpy2
        from rpy2.robjects import pandas2ri, Formula, conversion

        pandas2ri.activate()

        self.dds = None
        self.deseq_result = None
        self.comparison = None
        self.normalized_count_matrix = None
        self.gene_column = gene_column
        self.gene_id = count_matrix[self.gene_column]

        count_matrix = count_matrix.drop(gene_column, axis=1)

        logging.info(
            f'Number of columns in counts data {count_matrix.shape[1]} | '
            f'Number of rows in design matrix {design_matrix.shape[0]}'
        )

        # Load dataframe into R environment
        # Important: Change to r.data() if you use numpys and rpy2 latest versions
        self.count_matrix = conversion.py2rpy(count_matrix)

        # Assign columns to NULL
        self.count_matrix.names = rpy2.rinterface.NULL
        self.count_matrix = count_matrix
        self.design_matrix = conversion.py2rpy(design_matrix)
        self.design_formula = Formula(design_formula)

    def run(self, **kwargs):

        import rpy2
        from rpy2.robjects import pandas2ri, r, default_converter
        from rpy2.robjects.conversion import localconverter
        from rpy2.robjects.packages import importr

        pandas2ri.activate()
        deseq = importr('DESeq2')

        logging.info('imported DESeq2')

        to_dataframe = r('function(x) data.frame(x)')

        logging.info("running DESeq")

        self.dds = deseq.DESeqDataSetFromMatrix(
            countData=self.count_matrix,
            colData=self.design_matrix,
            design=self.design_formula
        )
        self.dds = deseq.DESeq(self.dds, **kwargs)
        # Previous script had "deseq.counts" instead
        # TODO check standard practice
        self.normalized_count_matrix = deseq.counts_DESeqDataSet(self.dds, normalized=True)

        self.comparison = deseq.resultsNames(self.dds)

        logging.info("Get DESeq2 results")

        self.deseq_result = deseq.results(self.dds, **kwargs)
        self.deseq_result = to_dataframe(self.deseq_result)
        self.deseq_result = rpy2.robjects.conversion.py2rpy(self.deseq_result)

        logging.info("Back to pandas dataframe")

        with localconverter(default_converter + pandas2ri.converter):
            results = rpy2.robjects.conversion.rpy2py(self.deseq_result)

        results[self.gene_column] = self.gene_id.values

        return results
