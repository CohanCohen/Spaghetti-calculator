import math
import numpy as np 
import matplotlib.pyplot as plt 
import os

def isbasic(string):
    for i in string:
        if (not i.isdigit() and i not in ["+", "-", "/", "*", "^", ".", "(", ")", "e"]):
            return False
    return True

def firstPr(string):
    k, j= 0, 0
    arr = []
    for i in string:
        if i == "(" and j == 0:k = 1
        if i == "(": j+=1
        if i == ")": j-=1
        if k == 1:
            arr.append(i)
        if j == 0 and k == 1:
            arr.pop(0);arr.pop(-1)
            return "".join(arr)
    return ""
def num(string):
    point = 0
    num_arr = [str(i) for i in range(0, 10)]
    num_arr.append("j")
    j = -1
    for i in string:
        j += 1
        if i not in num_arr:
            if i == ".": point += 1
            elif i == "-" and j != 0: return 0
            elif i == "-" : continue
            else: return 0
    if point > 1 : return 0
    return 1
def isspec(c):
    if 33<=ord(c)<=47 and ord(c) not in [40, 41]:
        return 1
    if ord(c) == 94:
        return 1
    if c == "_":
        return 1
    return 0
def isnumber(c):
    if 48<=ord(c)<=57:
        return 1
    return 0
def evalB(string):
    b = string.replace("^", "**")
    c = eval(b)
    if abs(c - round(c)) < 0.000000000001:
        return round(c)
    return c
def evalE(string, var, funcs):
    if isbasic(string):
        return evalB(string)
    else:
        n = string.find("asinh")
        if n != -1:
            pr = firstPr(string[n + 5:])
            expr = math.asinh(evalE(pr, var, funcs))
            nstr = string.replace("asinh("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("acosh")
        if n != -1:
            pr = firstPr(string[n + 5:])
            expr = math.acosh(evalE(pr, var, funcs))
            nstr = string.replace("acosh("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("atanh")
        if n != -1:
            pr = firstPr(string[n + 5:])
            expr = math.atanh(evalE(pr, var, funcs))
            nstr = string.replace("atanh("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("acoth")
        if n != -1:
            pr = firstPr(string[n + 5:])
            expr = math.atanh(1/evalE(pr, var, funcs))
            nstr = string.replace("acoth("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("sinh")
        if n != -1:
            pr = firstPr(string[n + 4:])
            expr = math.sinh(evalE(pr, var, funcs))
            nstr = string.replace("sinh("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("cosh")
        if n != -1:
            pr = firstPr(string[n + 4:])
            expr = math.cosh(evalE(pr, var, funcs))
            nstr = string.replace("cosh("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("tanh")
        if n != -1:
            pr = firstPr(string[n + 4:])
            expr = math.tanh(evalE(pr, var, funcs))
            nstr = string.replace("tanh("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("coth")
        if n != -1:
            pr = firstPr(string[n + 4:])
            expr = 1/math.tanh(evalE(pr, var, funcs))
            nstr = string.replace("coth("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("asin")
        if n != -1:
            pr = firstPr(string[n + 4:])
            expr = math.asin(evalE(pr, var, funcs))
            nstr = string.replace("asin("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("acos")
        if n != -1:
            pr = firstPr(string[n + 4:])
            expr = math.acos(evalE(pr, var, funcs))
            nstr = string.replace("acos("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("atan")
        if n != -1:
            pr = firstPr(string[n + 4:])
            expr = math.atan(evalE(pr, var, funcs))
            nstr = string.replace("atan("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("acot")
        if n != -1:
            pr = firstPr(string[n + 4:])
            expr = math.atan(1/evalE(pr, var))
            nstr = string.replace("acot("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("sin")
        if n != -1:
            pr = firstPr(string[n + 3:])
            expr = math.sin(evalE(pr, var, funcs))
            nstr = string.replace("sin("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("cos")
        if n != -1:
            pr = firstPr(string[n + 3:])
            expr = math.cos(evalE(pr, var, funcs))
            nstr = string.replace("cos("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("tan")
        if n != -1:
            pr = firstPr(string[n + 3:])
            expr = math.tan(evalE(pr, var, funcs))
            nstr = string.replace("tan("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("cot")
        if n != -1:
            pr = firstPr(string[n + 3:])
            expr = 1/math.tan(evalE(pr, var, funcs))
            nstr = string.replace("cot("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("log")
        if n != -1:
            pr = firstPr(string[n + 3:])
            expr = math.log10(evalE(pr, var, funcs))
            nstr = string.replace("log("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("ln")
        if n != -1:
            pr = firstPr(string[n + 2:])
            expr = math.log(evalE(pr, var, funcs))
            nstr = string.replace("ln("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("sqrt")
        if n != -1:
            pr = firstPr(string[n + 4:])
            expr = math.sqrt(evalE(pr, var, funcs))
            nstr = string.replace("sqrt("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("abs")
        if n != -1:
            pr = firstPr(string[n + 3:])
            expr = abs(evalE(pr, var, funcs))
            nstr = string.replace("abs("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("exp")
        if n != -1:
            pr = firstPr(string[n + 3:])
            expr = math.exp(evalE(pr, var, funcs))
            nstr = string.replace("exp("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("rad")
        if n != -1:
            pr = firstPr(string[n + 3:])
            expr = evalE(pr, var, funcs) * math.pi / 180
            nstr = string.replace("rad("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("floor")
        if n != -1:
            pr = firstPr(string[n + 5:])
            expr = math.floor(evalE(pr, var, funcs)) 
            nstr = string.replace("floor("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("ceil")
        if n != -1:
            pr = firstPr(string[n + 4:])
            expr = math.ceil(evalE(pr, var, funcs)) 
            nstr = string.replace("ceil("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        n = string.find("rnd")
        if n != -1:
            pr = firstPr(string[n + 3:])
            expr = round(evalE(pr, var, funcs)) 
            nstr = string.replace("rnd("+pr+")", "("+str(expr)+")")
            return evalE(nstr, var, funcs)
        for name, variable, st in funcs:
            n = string.find(name)
            if n != -1:
                pr = firstPr(string[n + len(name):])
                expr = evalE(str(evalE(st, [(variable, "("+str(evalE(pr, var, funcs))+")")], funcs)), var, funcs)
                nstr = string.replace(name+"("+pr+")", "("+str(expr)+")")
                return evalE(nstr, var, funcs)

        cond = False
        for i, j in var:
            if i in string:
                cond = True
            string = string.replace(i, "("+str(evalE(str(j), var, funcs))+")")  
        
        if cond : return evalE(string, var, funcs)
        return 0
    
def ins_star(string):
    new_string = string[:]
    k = len(new_string) - 1
    i = 0
    while i < k:
        if isnumber(new_string[i]) + isnumber(new_string[i+1]) == 1:
            if new_string[i+1] != ")" and new_string[i] != "(" :
                if isspec(new_string[i]) + isspec(new_string[i+1]) == 0:
                    new_string = new_string[:i+1] + "*" + new_string[i+1:]
                    k+=1
        if new_string[i] == ")" and 97<=ord(new_string[i+1])<=122:
            new_string = new_string[:i+1] + "*" + new_string[i+1:]
            k+=1
        i += 1
    return new_string
def norm_pl(string, var, funcs):# mononominals_normed = [[pow, coef], [pow, coeff], ...]
    nstring = ins_star(string)
    mononomials = nstring.split("+")
    mononomials_normed = []
    for mn in mononomials:
        ml = evalE(mn, [(var, "1")], funcs)
        pow = round(math.log(evalE("("+mn+")/"+str(ml), [(var, str(math.e))], funcs)))
        mononomials_normed.append([pow, ml])
    return mononomials_normed


def nDiff(function, n, dx=0.00001):
    if n == 1:
        return lambda x : (function(x+dx) - function(x))/dx
    else:
        f = lambda x : (function(x+dx) - function(x))/dx
        return nDiff(f, n - 1)
    
def nDiffStr(string, var, funcs,n, dx=0.00001): #var = "x"
    f = lambda x : evalE(string, [(var, x)], funcs)
    return nDiff(f, n, dx=dx)

def nInt(function, a, b, dx=0.0000001):
    s = 0
    i = a
    while i < b:
        s += function(i) * dx
        i += dx
    return s

def nIntStr(string, var, funcs, a, b, dx=0.00001): # e.g var = "x"
    f = lambda x : evalE(string, [(var, x)], funcs)
    return nInt(f, a, b, dx=dx)

def plotfunction(function, domain):
    X = []
    Y = []
    i = domain[0]
    while i < domain[1]:
        X.append(i)
        Y.append(function(i))
        i += 0.01
    min_x, max_x = min(X), max(X)
    min_y, max_y = min(Y), max(Y)
    plt.autoscale(False)
    ax = plt.gca()
    ax.set_xlim([min_x, max_x])
    ax.set_ylim([min_y, max_y])
    plt.arrow((min_x + max_x)/2, min_y, 0, max_y-min_y)
    plt.arrow(min_x, 0, max_x-min_x, 0)
    plt.plot(X, Y)
    plt.show()

def plotfunctionStr(string, var, vars, funcs, domain):
    f = lambda x : evalE(string, [(var, x)]+vars, funcs)
    plotfunction(f, domain)

class polynomial:
    def __init__(self, arr):
        narr = {}
        for p, c in arr:
            if p not in narr.keys():
                narr.update({p : c})
            else:
                narr[p] += c 
        
        self.coeff_arr = sorted(list(narr.items()), key=lambda x:x[0], reverse=True)
        self.deg = max([i for i,j in arr])
    
    def __add__(self, other):
        arr = []
        for p1, c1 in self.coeff_arr:
            for p2, c2 in other.coeff_arr:
                if p1 == p2:
                    arr.append([p1, c1+c2])
        return polynomial(arr)
    
    def __sub__(self, other):
        arr = []
        for p1, c1 in self.coeff_arr:
            for p2, c2 in other.coeff_arr:
                if p1 == p2:
                    arr.append([p1, c1-c2])
        return polynomial(arr)
    
    def __mul__(self, other):
        arr = []
        for p1, c1 in self.coeff_arr:
            for p2, c2 in other.coeff_arr:
                arr.append([p1 + p2, c1*c2])
        return polynomial(arr)
    
    def __call__(self, x):
        s = 0
        for p, c in self.coeff_arr:
            s += c * x**p
        return s
    
    def __str__(self):
        string = ""
        for p, c in self.coeff_arr:
            if c != 1:
                if p > 1 : string += str(c) + "x^" + str(p)+"+" 
                if p == 1 : string += str(c) + "x+"
                if p == 0 : string += str(c)
            else:
                if p > 1 : string += "x^" + str(p)+"+" 
                if p == 1 : string += "x+"
                if p == 0 : string += str(c)
        return string[:-1]

    def diff(self):
        new_arr = []
        for p, c in self.coeff_arr:
            if p > 0:
                new_arr.append([p - 1, c * p])
        
        return polynomial(new_arr)
    
    def root(self):
        if self.deg == 0:
            return 0
        if self.deg == 1:
            for p, c in self.coeff_arr:
                if p == 1:
                    a = c 
                if p == 0:
                    b = c
            return -b/a 
        if self.deg == 2:
            for p, coeff in self.coeff_arr:
                if p == 2:
                    a = coeff 
                if p == 1:
                    b = coeff
                if p == 0:
                    c = coeff
            if b**2-4*a*c < 0:
                return None
            return (-b + math.sqrt(b**2-4*a*c)) / (2*a)
        
        if self.deg == 3:
            for p, coeff in self.coeff_arr:
                if p == 3:
                    a = coeff 
                if p == 2:
                    b = coeff
                if p == 1:
                    c = coeff
                if p == 0:
                    d = coeff
            f0 = (-b**3)/(27*a**3) + b*c/(6*a**2)-d/(2*a)
            f1 = f0**2+(c/(3*a) - (b**2)/(9*a**2))**3
            if f1<0:
                return None 
            return (f0+math.sqrt(f1))**(1/3)+(f0-math.sqrt(f1))**(1/3)-b/(3*a)
                  
        return None       

cond = True
pm = False
variables = [("pi", math.pi), ("E", math.e)]
functions = []
polynomials = {}
while cond : 
    comm = input(">>> ").lower().split(" ")
    if comm[0] == "quit":
        cond = False
        break
    elif comm[0] == "eval":
        print("--> ", evalE(ins_star(comm[1]), variables[:], functions[:]))
    elif comm[0] == "var":
        variables.append([comm[1], comm[2]])
    elif comm[0] == "let":
        c = comm[1].split("=")
        v = firstPr(c[0])
        n = c[0][:c[0].find("(")]
        s = ins_star(c[1])
        functions.append([n, v, s])
    elif comm[0] == "drvt" : # drvt x f(x) 5
        print("--> ", nDiffStr(ins_star(comm[2]), comm[1], functions[:], 1)(evalE(ins_star(comm[3]), variables[:], functions[:])))
    
    elif comm[0] == "graph":# graph 2f(x) -1 1
        print("->notate variable via x.")
        plotfunctionStr(ins_star(comm[1]), "x", variables[:], functions[:], [evalE(comm[2], variables[:], functions[:]), evalE(comm[3], variables[:], functions[:])])
    
    elif comm[0] == "int":
        print("--> ", nIntStr(ins_star(comm[2]), comm[1], functions[:], evalE(ins_star(comm[3]), variables[:], functions[:]), evalE(ins_star(comm[4]), variables[:], functions[:])))
    
    elif comm[0] == "help":
        with open("doc.txt", "r") as f:
            print("->\n")
            print(f.read())
            f.close()
    elif comm[0] == "setpl":
        c = comm[1].split("=")
        v = firstPr(c[0])
        n = c[0][:c[0].find("(")]
        s = ins_star(c[1])
        polynomials.update({n: polynomial(norm_pl(c[1], v, functions[:]))})
        functions.append([n, v, s])
    elif comm[0] == "-pm":
        pm = True 
    
    elif pm:
        if comm[0] == "mul":
            p3 = polynomial([[0, 1]])
            j = 0
            for i in comm:
                if j <= 1:
                    continue
                p3 *= polynomials[i]
            polynomials.update({comm[1] : p3})
            functions.append([comm[1], "x", str(p3)])
            print("--> ", p3)
        if comm[0] == "add":
            p3 = polynomial([[0, 0]])
            j = 0
            for i in comm:
                if j <= 1:
                    continue
                p3 += polynomials[i]
            polynomials.update({comm[1] : p3})
            functions.append([comm[1], "x", str(p3)])
            print("--> ", p3)
        if comm[0] == "rt":
            print("-->", polynomials[comm[1]].root())
    else:
        print("->Invalid command type help for help.\n")

quit()
    