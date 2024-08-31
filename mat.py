def matmul(a, b):
    result = [[0 for i in range(len(a))] for j in range(len(b[0]))]
    for i in range(len(a)):
        for j in range(len(b[0])):
            for k in range(len(b)):
                result[i][j] += a[i][k] * b[k][j]
    return result
x = [[1,2],[4,5],[7,8]]
y = [[1,2,3],[4,5,6]]
print(matmul(x,y))
