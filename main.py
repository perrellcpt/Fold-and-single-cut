import numpy as np
import scipy
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import svgwrite
import sys
import time
import copy


def vprint(v):  # given a list, print all elements in the list
    for x in range(0, len(v)):
        print(v[x])

def ccw(A, B, C):
    return (C[1] - A[1]) * (B[0] - A[0]) > (B[1] - A[1]) * (C[0] - A[0])


# Return true if line segments AB and CD intersect
def intersect(A, B, C, D):
    return ccw(A, C, D) != ccw(B, C, D) and ccw(A, B, C) != ccw(A, B, D)


def plot(vert):  # given list of ordered pairs of verticies, plot lines
    X = []
    Y = []
    for x in range(0, len(vert)):
        X.append(vert[x][0])
        Y.append(vert[x][1])

    X.append(vert[0][0])
    Y.append(vert[0][0])
    x1 = np.array(X)
    y1 = np.array(Y)
    for n in range(1, len(X) - 1):
        for m in range(n + 1, len(X)):
            if intersect([X[n - 1], Y[n - 1]], [X[n], Y[n]], [X[m - 1], Y[m - 1]], [X[m], Y[m]]) == True:
                print("invalid polygon, lines intersect")
                exit(2)
    plt.plot(x1, y1)


def nonincident(v, E):  # given a vertex and the list of edges, remove the incident edges
    nonI = copy.deepcopy(E)
    x=len(nonI)-1
    while x>=0:
        if v==E[x][0] or v == E[x][1]:
            nonI.pop(x)
        x-=1
    return nonI


def diskpack1(V, E):
    # given a list of ordered vertices and a list of edges, define the radius of a circle that is halfway between the
    # vertex and the nearest non-incident edge for every vertex in the list and return a list of vertex-radius pairs
    rad = []
    for v in range(0, len(V)):
        n=-1
        #print("the vertex is:")
        #vprint(V[v])
        #print("calling nonincident:\n")
        NI = nonincident(V[v], E)
        #print("the list of edges without the vertex: \n")
        #vprint(NI)
        vert = np.array(V[v])
        for e in range(0, len(NI)):

            norm = np.array([NI[e][0][1] - NI[e][1][1], -(NI[e][0][0] - NI[e][1][ 0])])  # normal to the vector from
            # the point NI[e][1] to the point NI[e][2]
            A = np.array(NI[e][0])
            B = np.array(NI[e][1]) #vectors of endpoints of the edge
            t=(np.dot((vert-A), (B-A))) / ((np.linalg.norm(B-A))**2) #calculate the projection parameter
            if t<0:
                r=abs(np.linalg.norm(vert-A)) #the closest point on the edge to the vertex is A
            elif t>1:
                r=abs(np.linalg.norm(vert-B)) #the closest point on the edge to the vertex is B
            else: #otherwise, assume the edge to be an infinite line and use vector projection
                p = np.array([V[v][0] - NI[e][0][0], V[v][1] - NI[e][0][1]])  # vector from one end of the edge to the vertex
                r= abs(np.dot(p, norm) / np.linalg.norm(norm)) #vector projection
            #print("the vertex is " + str(r) + " distance away from the edge ")
            #vprint(NI[e])
            if r<n or n==-1: #want the distance of the closest edge to the vertex
                n=r
        rad.append([V[v], n/2])
    return rad

def diskpack2(E,V,R, cov, ucov): #given a list of edges, vertices and the corresponding list of radii of the circles placed on each vertex, split each edge on the point
    # where the edge of the circle intersects with the edge
    for x in range(0,len(E)):
        ucov.append([])
        for y in range(0,2):
            print("edge: \n")
            vprint(E[x])
            v = np.array([-E[x][y][0]+E[x][(y+1)%2][0], -E[x][y][1]+E[x][(y+1)%2][1]]) #vector from one endpoint to another
            print("vector: \n")
            vprint(v)
            v = (v/np.linalg.norm(v))
            for z in range(0, len(R)):
                if R[z][0] == E[x][y]:
                    v = v * R[z][1]
                    print("ok\n")
            v= v + np.array(E[x][y])
            print("modified vector: \n ")
            vprint(v)
            cov.append([E[x][y], [float(v[0]), float(v[1])]])
            ucov[x].append([float(v[0]),float(v[1])])

def diskpack3(ucov, R): # take the uncovered edges and cover them with diameter disks, splitting if disks overlap
    crowded=[]
    midp=[]
    leng=[]
    skip=0
    done=[]
    while 1:
        for x in range(0, len(ucov)):
            midp.append([(ucov[x][0][0]+ucov[x][1][0])/2, (ucov[x][0][1]+ucov[x][1][1])/2]) #midpoint of the edge
            leng.append(abs(np.linalg.norm(np.array(ucov[x][0]) - np.array(ucov[x][1])))) # length of the edge
        for x in range(0,len(midp)):
            for y in range(0, len(midp)):
                if x == y:
                    continue
                dist=abs(np.linalg.norm(np.array(midp[x]) - np.array(midp[y]))) #calculate the distance between the two points
                if dist < (leng[x] + leng[y]):
                    crowded.append([ucov[x][0], midp[x]]) #if the diameter disks overlap, split the edge and put it in crowded
                    crowded.append([ucov[x][1], midp[x]])
                    skip=1
                    break
            if skip==1:
                skip=0
                continue
            for y in range(0, len(R)):
                dist=abs(np.linalg.norm(np.array(midp[x]) - np.array(R[y][0])))
                if dist < (leng[x] + R[y][1]):
                    crowded.append(x)
                    skip = 1
                    break
            if skip==0:
                done.append(ucov[x])

        if len(crowded)==0:
            break
        ucov.clear()
        ucov=copy.deepcopy(crowded)
        crowded.clear()
        midp.clear()
        leng.clear()
    for x in range(0, len(done)):
        c = pat.Circle((done[x][0][0]+done[x][1][0])/2, (done[x][0][1]+done[x][1][1])/2,
                       radius=abs(np.linalg.norm(np.array(done[x][0]) - np.array(done[x][1]))), color='yellow')
        ax.add_patch(c)


num = int(input("Enter the number of points in the polygon: "))
if num < 2:
    print("not enough vertices were given")
    sys.exit(1)
print("Enter the points of the polygon in order in the form x,y:")
vertices = []
for x in range(0, num):
    coords = input().split(",")
    coords[0] = float(coords[0])
    coords[1] = float(coords[1])
    if (coords[0] > 1 or coords[0] < 0) or (coords[1] > 1 or coords[1] < 0):
        print("invalid locations for vertices")
        exit(1)
    vertices.append(coords)

vprint(vertices)
corners = [[0, 0], [0, 1], [1, 1], [1, 0]]
plot(corners)
plot(vertices)


edges = []

for x in range(1, len(vertices)):
    edges.append([vertices[x - 1], vertices[x]])

edges.append([vertices[0],vertices[len(vertices)-1]])
for x in range(1, len(corners)):
    edges.append([corners[x - 1], corners[x]])

edges.append([corners[0],corners[len(corners)-1]])
#vprint(edges)

vertices.extend(corners) #all the original vertices in PR
R=diskpack1(vertices, edges)
ax=plt.subplot()
# c=pat.Circle((0.5,0.5), radius=0.3, color='blue')
# ax.add_patch(c)

for x in range(0, len(vertices)):
    c = pat.Circle((float(vertices[x][0]), float(vertices[x][1])), radius=R[x][1], color='blue')
    ax.add_patch(c)

UcovE=[] #uncovered edges
covE=[] #covered edges
diskpack2(edges,vertices,R, covE, UcovE)
vprint(covE)
print("goy \n")
vprint(UcovE)
plt.show()
