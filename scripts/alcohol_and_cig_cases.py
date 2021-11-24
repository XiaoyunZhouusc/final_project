alcohol_cases = None
with open('../static/alcohol_cases.txt', 'r') as f:
    alcohol_cases = f.readlines()

smoke_cases = None
with open('../static/smoke_cases.txt', 'r') as f:
    smoke_cases = f.readlines()

with open('../static/alcohol_and_cig_cases.txt', 'w') as f:
    for i in smoke_cases:
        if i in alcohol_cases:
            f.write(i)
