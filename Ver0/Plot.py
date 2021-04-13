from matplotlib import pyplot

with open("O.txt", "r") as O:
    Lines = O.readlines()

Squares = [list(map(float, Line.strip().split())) for Line in Lines]

pyplot.figure(figsize = (8, 8))
pyplot.xlim(-0.25, 1.25); pyplot.ylim(-0.25, 1.25)
for Sq in Squares:
    x, y, s = Sq
    pyplot.plot([x-s/2, x+s/2, x+s/2, x-s/2, x-s/2], [y-s/2, y-s/2, y+s/2, y+s/2, y-s/2], '-k')
pyplot.show()
pyplot.close()