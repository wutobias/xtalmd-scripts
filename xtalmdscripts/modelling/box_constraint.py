# constraint for trust-constr minimization to avoid openFF Periodic Boundary Conditions
def constraint1(x):
    box_vector = x[-6:]
    value1 = box_vector[0] - 2 * abs(box_vector[1])
    return value1

def constraint2(x):
    box_vector = x[-6:]
    value2 = box_vector[0] - 2 * abs(box_vector[3])
    return value2

def constraint3(x):
    box_vector = x[-6:]
    value3 = box_vector[2] - 2 * abs(box_vector[4])
    return value3

