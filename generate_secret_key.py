from django.core.management import utils

if __name__ == '__main__':
    print(utils.get_random_secret_key())
