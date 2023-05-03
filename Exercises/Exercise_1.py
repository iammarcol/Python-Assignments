# PROBLEM_1

def get_sphere_volume(radius):
    "Calculate the volume of a sphere"
    return float((4*3.142*radius**3)/3)

# PROBLEM_2a

def recursive_factorial(n):
    "Given the input num n the function calculates its factorial"
    if n==1:
        return 1
    else:
        return (n*recursive_factorial(n-1))


# PROBLEM_2b

def factorial(n):
    fact=1
    i=1
    if n>=0:
        while 0<i<=n:
            fact=fact*i
            i+=1
    return fact

# PROBLEM_3a

def recursive_count_up(n,odd=''):
    "If 'True' is not specified odd is set to false"
    if odd==True:
        if n>=0:
            recursive_count_up(n-1,odd)
            if n%2!=0:
                print(n)
    else:
        if n>=0:
            recursive_count_up(n-1,odd)
            print(n)

# PROBLEM_3b

def count_up(n, odd=''):
    "If 'True' is not specified odd is set to false"
    for i in range(0, n+1):
        if odd==True:
            if i%2!=0:
             print(i)
        else:
             print(i)


# PROBLEM_4

def get_final_price(price, discount_percentage=10):
    """Return the final price after applying the discount percentage"""
    return (price-(price*discount_percentage/100))
# if only one argument is given, discount percentage is set to 10 by default





