#!/usr/bin/env python3

import pandas as pd
import sys, decimal, re
from tkinter import Tk

def getmeasurement(x_input, dx_input) -> str:
    """
        Get shorthand notation for a measurement with uncertainty: 1.213(3).
    """
    ctx = decimal.Context()
    ctx.rounding = decimal.ROUND_HALF_UP
    x = decimal.Decimal(x_input)
    dx = decimal.Decimal(dx_input)
    if dx == decimal.Decimal('0'):
        return repr(x_input)

    #dx = ctx.create_decimal(repr(dx_input)
    #dx = ctx.create_decimal(dx_input)
    #ctx.prec = 20
    #std_str = format(dx, 'f')
    #print(std_str)

    # obvious check
    if dx < decimal.Decimal('0'):
        return [str(x), str(dx)]

    # value to add at the end for the precision
    add = 1

    # distinguish between intervals of uncertainty
    if decimal.Decimal('1') <= dx and dx < decimal.Decimal('19.5'):
        digits = re.findall(r"[0-9]", str(dx))
        if dx >= decimal.Decimal('10'):
            add = 0
    elif dx < decimal.Decimal('1'):
        if re.search(r"-", str(dx)) is not None:
            #mantissa, exponent = str(dx).replace("E", "").split('-')
            print("Exponent uncertainty E-7 too small.")
            return
    else:
        print("Found uncertainty > 19. To avoid confusion, use scientific notation.")
        return

    # check if there's the dot for decimal places
    dot = re.search(r"\.", str(dx))
    if dot is not None:
        string_dx = str(dx).split('.')[1].lstrip('0')
        if dx >= decimal.Decimal('1'):
            string_dx = str(dx).replace(".", "")
        digits = re.findall(r"[0-9]", string_dx)
    
    # uncertainty to first digit
    dx_str = digits[0]

    # check to round uncertainty
    first_digit_reg = re.search(r"[1-9]", str(dx))
    first_digit = first_digit_reg.group(0)
    if first_digit != '1':
        # check if the following digit is > 4, then round
        if len(digits) >= 2 and digits[1] in [f'{i}' for i in range(5,10)]:
            dx_str = str(int(digits[0])+1)
    else:
        # check if first digit of uncertainty is one, if so add another digit; if there's none other, add "0"
        add += 1
        last_digit = "0"
        if len(digits) >= 2:
            last_digit = digits[1]
            if len(digits) >= 3 and digits[2] in [f'{i}' for i in range(5,10)]:
                if digits[1] != '9':
                    last_digit = str(int(digits[1])+1)
                else:
                    dx_str = ''
                    last_digit = '2'
                    add -= 1
            dx_str += last_digit

    
    # precision to round at
    prec = first_digit_reg.start()-2+add

    # edge cases
    if prec == 0:
        prec = 1
    if prec < 0:
        prec = 0

    #return [str(x), dx_input, f"{x:.{prec}f}", dx_str, first_digit_reg.start(), add, prec]
    #return [f"{x:.{prec}f}", dx_str]
    return f"{x:.{prec}f}({dx_str})"

def csvtotables(file_name, dati):
    # read x, y and y error from file using run number
    col_list = [f"x {dati}", f"xerr {dati}", f"y {dati}", f"yerr {dati}"]
    df = pd.read_csv(file_name, usecols=col_list, sep=",", decimal=".")#, dtype=str)#, converters={i: str for i in range(0,16)})

    # manipulate data into array
    df = df.dropna()
    df = df.sort_values([df.columns[0], df.columns[2]]).reset_index()
    xdata = df[f"x {dati}"]
    ydata = df[f"y {dati}"]
    xerr = df[f"xerr {dati}"]
    yerr = df[f"yerr {dati}"]

    #file = open("data.txt", "w")
    #for i in range(len(xdata)):
    #    file.write(repr(xdata[i]) + " " + repr(ydata[i]) + " " + repr(xerr[i]) + " " + repr(yerr[i]) + "\n")
    #file.close
    clipboard = Tk()
    clipboard.withdraw()
    clipboard.clipboard_clear()
    for i in range(0,len(xdata)):
        s_x = getmeasurement(xdata[i], xerr[i])
        s_y = getmeasurement(ydata[i], yerr[i])
        clipboard.clipboard_append("\\SI{"+ s_x +"}{}" + " & " + "\\SI{"+ s_y +"}{}" + " \\\\" + "\n")
    clipboard.update()

if __name__ == '__main__':
    if (len(sys.argv) == 4):
        print(getmeasurement(*sys.argv[1:3]))
    else:
        csvtotables(*sys.argv[1:])