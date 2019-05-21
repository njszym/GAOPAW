import sys


## Just placeholders for now
diff = 0
tol = 1

## Need to evaluate difference in valence energy from AtomPAW at each step

while diff < tol:
    filename = sys.argv[-1]
    with open(filename) as f:
        lines = f.readlines()
    log_line = (lines[1]).split()
    index = 1
    for elem in log_line:
        if elem == 'loggrid':
            num_pts = float(log_line[index])
            break
        else:
            index += 1
    print (num_pts)
    num_pts -= 10
    log_line[index] = str(num_pts)
    new_line = ''
    for elem in log_line:
        var = elem+' '
        new_line += var
    lines[1] = new_line+'\n'
    with open(filename,'w') as f:
        for line in lines:
            f.write(line)
    diff += 0.1
