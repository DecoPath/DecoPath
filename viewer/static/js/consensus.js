/*! Consensus table styling */

function addConsensusColour(consensus_value) {
    if (consensus_value == 0) {
        return 'Discordant';
    } else if (consensus_value == 1) {
        return 'NoMappings';
    } else return 'Concordant';
}

function getConsensus(table_consensus, consensus_value) {
    var colour = addConsensusColour(consensus_value);
    if (colour == 'NoMappings') {
        return '<a class="btn btn-primary ' + colour + '">No mappings</a>';
    } else return '<a class="btn btn-primary ' + colour + '">' + colour + '</a>';
}

// Add DC-DB consensus comparison checker
function addDCConsensus(dc_val) {
    if (dc_val == 0) {
        return "discordant"
    } else if (dc_val == 2) {
        return "concordant"
    } else return "no-mappings"
}

// Add colour to score when over-expressed (pink) and under-expressed (blue) for ORA
function addColourOra(qval) {
    if (qval > qval_threshold) {
        return 'not-significant';
    } else return 'significant-ora';
}

// Add colouring to score when over-expressed (pink), under-expressed (blue) and not-significant (gray) for GSEA
function addColourGsea(score, qval) {
    if (qval > default_qval) {
        return 'not-significant';
    } else if (score <= -2) {
        return 'high-under-expressed';
    } else if (score < -1 && score > -2) {
        return 'under-expressed';
    } else if (score >= -1 && score < 0) {
        return 'low-under-expressed';
    } else if (score == 0) {
        return 'no-change';
    } else if (score > 0 && score <= 1) {
        return 'low-over-expressed';
    } else if (score > 1 && score < 2) {
        return 'over-expressed';
    } else if (score >= 2) {
        return 'high-over-expressed';
    }
}

// Add link to score ORA
function addLinkOra(identifier, qval, dc_val) {

    var colour = addColourOra(qval);
    var dc_colour = addDCConsensus(dc_val)

    if (identifier.startsWith("hsa")) {
        return '<a href="https://www.genome.jp/dbget-bin/www_bget?pathway+' + identifier + '" class="btn btn-primary ' + colour + '" target="_blank" >' + qval + '</a>';
    } else if (identifier.startsWith("PW")) {
        return '<a href="https://pathbank.org/pathwhiz/pathways/' + identifier + '" class="btn btn-primary ' + colour + '" target="_blank">' + qval + '</a>';
    } else if (identifier.startsWith('R-HSA')) {
        return '<a href="https://reactome.org/PathwayBrowser/#/' + identifier + '" class="btn btn-primary ' + colour + '" target="_blank">' + qval + '</a>';
    } else if (identifier.startsWith("WP")) {
        return '<a href="https://www.wikipathways.org/index.php/Pathway:' + identifier + '"class="btn btn-primary ' + colour + '" target="_blank">' + qval + '</a>';
    } else if (identifier.startsWith("DC")) {
        return '<p><span style="float:left;" class="dot ' + dc_colour + '"></span><a href="../dc_genesets/' + identifier + '" class="btn btn-primary ' + colour + '" target="_blank">&nbsp;' + qval + '</a></p>';
    } else if (identifier == 'NA') {
        return '<p class="no-mapping" style="text-align: left">No mapping</p>';
    } else return '<a class="btn btn-primary ' + colour + '" target="_blank">' + qval + '</a>';
}


// Add link to score GSEA
function addLinkGsea(identifier, score, qval, dc_val) {

    var colour = addColourGsea(score, qval);
    var dc_colour = addDCConsensus(dc_val)

    if (identifier.startsWith("hsa")) {
        return '<a href="https://www.genome.jp/dbget-bin/www_bget?pathway+' + identifier + '" class="btn btn-primary ' + colour + '" target="_blank" >' + score + '</a>';
    } else if (identifier.startsWith("PW")) {
        return '<a href="https://pathbank.org/pathwhiz/pathways/' + identifier + '" class="btn btn-primary ' + colour + '" target="_blank">' + score + '</a>';
    } else if (identifier.startsWith("R-HSA")) {
        return '<a href="https://reactome.org/PathwayBrowser/#/' + identifier + '" class="btn btn-primary ' + colour + '" target="_blank">' + score + '</a>';
    } else if (identifier.startsWith("WP")) {
        return '<a href="https://www.wikipathways.org/index.php/Pathway:' + identifier + '"class="btn btn-primary ' + colour + '" target="_blank">' + score + '</a>';
    } else if (identifier.startsWith("DC")) {
        return '<p><span style="float: left" class="dot ' + dc_colour + '"></span><a href="../dc_genesets/' + identifier + '" class="btn btn-primary ' + colour + '" target="_blank">&nbsp;' + score + '</a></p>';
    } else if (identifier == "NA") {
        return '<p class="no-mapping" style="text-align: left">No mapping</p>';
    } else return '<a class="btn btn-primary ' + colour + '" target="_blank">' + score + '</a>';
}
