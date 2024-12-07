from datetime import datetime
from math import floor
from time import sleep


def plotBar(x):
    bar = f'|{"█" * floor(x * 50):<50}|'
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> {bar}.', end = '\r' if x < 1 else '\n', flush = True)
    return None

if __name__ == '__main__':
    for i in range(1, 101):
        sleep(0.5)
        plotBar(i / 100)
