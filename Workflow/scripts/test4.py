ok="3S6M1I6M1I7M1I10M1I24M2D2M1I3M1D5M1D32M2D7M1I13M1D4M1D3M1D5M1D14M1D5M1D24M1D11M1I16M1D2M2D11M2D5M1D30M1D20M1I14M1D2M1D2M1D3M1I19M2D4M1D6M1I18M1D4M1D4M1D5M1D17M1D4M1I17M1D3M1I4M2D29M1D1M1D9M1D10M1D11M1D10M1D22M1I5M1D14M1D60M1I26M1I16M1D21M1I1M1I33M1D7M1D28M1I8M1I13M1D20M1D10M1I20M1D8M1I4M1D6M1D1M2I11M1D13M1I5M1I6M1I30M1D14M1D34M1D1M1D2M1I44M1I13M1I9M1D32M1D6M2I38M1D4M1D22M1D33M1D7M1D5M1I7M1D9M1D35M1D3M1D20M1I32M1I8M2D1M2D4M2I5M1I41M1D1M1D26M1D33M1D45M4S"
rep=''
seq=''
for i in range(len(ok)):
    if ok[i].isnumeric():
        rep+=ok[i]
        print(rep)
    else:
        if ok[i]=='M':
            car='='
            seq+=int(rep)*car
            rep=''
        if ok[i]=='S':
            car='.'
            seq+=int(rep)*car
            rep=''
        if ok[i]=='I':
            car='*'
            seq+=int(rep)*car
            rep=''
        if ok[i]=='D':
            car='_'
            seq+=int(rep)*car
            rep=''
        if ok[i]=='N':
            car='?'
            seq+=int(rep)*car
            rep=''
        if ok[i]=='H':
            car=':'
            seq+=int(rep)*car
            rep=''
        if ok[i]=='P':
            car='!'
            seq+=int(rep)*car
            rep=''
print(seq)