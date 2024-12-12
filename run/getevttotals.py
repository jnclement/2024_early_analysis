njt = 0
eff = 0

with open("runevt.txt","r") as thefile:
    for line in thefile:
        line = line.strip().split()
        if len(line) > 1:
            eff += float(line[5])
            njt += int(line[4])

print(eff)
print(njt)
print((eff/(21*(10**-3)))/(10**12))
