def vprint(v):
    for x in range(0, len(v)):
        print(v[x])

num=int(input("Enter the number of points in the polygon: "))
print("Enter the points of the polygon in order in the form x,y:")
vertices=[]
for x in range(0,num):
    coords=input().split(",")
    coords[0]=int(coords[0])
    coords[1]=int(coords[1])
    vertices.append(coords)

vprint(vertices)
