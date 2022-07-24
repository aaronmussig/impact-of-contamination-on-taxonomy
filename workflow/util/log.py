from datetime import datetime


def log(msg: str, title=False):
    if title:
        print('\n' + '-' * 80)
    print(f'[{datetime.now().strftime("%d/%m/%Y %H:%M:%S")}] - {msg}')
    if title:
        print('-' * 80)
