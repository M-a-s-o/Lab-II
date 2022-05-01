import decimal, re, numpy as np
from sys import argv, exit

def auto_std(data):
    mean = data[0]
    std = data[1]
    ctx = decimal.Context()
    d1 = ctx.create_decimal(repr(std))
    ctx.prec = 20
    std_str = format(d1, 'f')

    n_decimals = 0
    match = re.search("[1-9]", std_str)
    if(match.group() == "1" or match.group() == "2"):
        mean = np.round(mean , match.start())
        std = np.round(std , match.start())
    else:
        mean = np.round(mean , match.start() - 1)
        std = np.round(std , match.start() - 1)
  
    return np.array([mean , std])


if __name__ == '__main__':
    if (len(argv) != 3):
        print("Required two arguments.")
        exit(1)
    arr = [float(argv[i]) for i in range(1,3)]
    print(auto_std(arr))