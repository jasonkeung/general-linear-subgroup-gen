# matrix mult function
def matMult(m1, m2, mod):
    tl = m1[0] * m2[0] + m1[1] * m2[2]
    tr = m1[0] * m2[1] + m1[1] * m2[3]
    bl = m1[2] * m2[0] + m1[3] * m2[2]
    br = m1[2] * m2[1] + m1[3] * m2[3]

    return (tl % mod, tr % mod, bl % mod, br % mod)

# inverse function mod p
def invMod(num, mod):
    for i in range(mod):
        if (num * i) % mod == 1:
            return i

    return -1 # this should never happen if mod is prime and num is nonzero

# matrix print function
def printMat(mat):
    print('|' + str(mat[0]) + ' ' + str(mat[1]) + '|')
    print('|' + str(mat[2]) + ' ' + str(mat[3]) + '|')
    print()

# matrix determinant function
def matDet(mat):
    return mat[0] * mat[3] - mat[1] * mat[2]

print('Welcome to Jason\'s general linear group subgroup generator!')
mod = int(input('Enter a prime number p to use GL(Z_p) as your group: '))

print('Matrix multiplication mod ' + str(mod) + ' is the operation for this general linear group.')
stack = []
subgroup = set()
subgroup.add((1, 0, 0, 1)) # ensure identity element is in subgroup

print('Enter an invertible 2x2 matrix from left-right, top-down to generate its subgroup: ')

tl = int(input('   Top  left: '))
tr = int(input('  Top  right: '))
bl = int(input('Bottom  left: '))
br = int(input('Bottom right: '))

if matDet([tl, tr, bl, br]) != 0:
    stack.append((tl % mod, tr % mod, bl % mod, br % mod))


while stack: # stack will hold elements to be processed
    a = stack.pop() # element to process
    
    if matDet(a) != 0:
    
        if a not in subgroup:
            subgroup.add(a)

        # add inverse of a to stack:
        detInv = invMod(matDet(a), mod) # 1/determinant in mod p
        matSwapNeg = (a[3], -a[1], -a[2], a[0]) # swap a and d, negate b and c
        matInv = (detInv * matSwapNeg[0] % mod, detInv * matSwapNeg[1] % mod, detInv * matSwapNeg[2] % mod, detInv * matSwapNeg[3] % mod) # distribute factor in
        if matInv not in stack and matInv not in subgroup:
            stack.append(matInv)



        # add a * b and b * a for every b already in subgroup
        for b in subgroup:
            prod1 = matMult(a, b, mod) # a * b
            if prod1 not in stack and prod1 not in subgroup:
                stack.append(prod1)

            prod2 = matMult(b, a, mod) # b * a
            if prod2 not in stack and prod2 not in subgroup:
                stack.append(prod2)

print('Subgroup generated with ' + str(len(subgroup)) + ' elements:')
for mat in subgroup:
    printMat(mat)


    
