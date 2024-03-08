# MCS 595 Thu 23 Jan 2013 : drawtrop.py
"""
This script draws a plane tropical curve.
0. The main function prompts the user for a degree d.
   Then all exponents for monomials up to the given degree d are computed. 
   A random number generator makes the corresponding list of random
   coefficients in [-1, +1].
1. The drawing of the curve happens in two stages.
   First we compute the nodes, that are those points where the minimum
   is attained (at least) three times.  
   The script generates all triplets of indices to the monomials.
   For each triplet, the intersection point is computed and then
   those points are retained as nodes after testing whether the
   minimum is indeed attained three times at the point.
2. In the second stage the lines connecting the nodes are mades.
   Those lines are those places where the minimum is attained twice.
   We distinguish between edges between the nodes and half rays.
   An edge connects two nodes.  To decide whether two nodes are connected, 
   we consider the midpoint and check whether the minimum is attained
   twice at that midpoint.
   At each node we then draw the half rays, along the directions of
   each pair of lines that intersect at the node.
"""

import matplotlib.pyplot as plt

TOL = 1.0e-4
OUTMAX = 100

def exponents(deg):
    """
    Generates a list of 2-tuples with the exponents
    of all monomials up to degree deg.
    """
    deg = int(deg)
    xrd = range(deg+1)
    yrd = range(deg+1)
    result = []
    for degx in xrd:
        for degy in yrd:
            if(degx + degy <= deg):
                result.append((degx, degy))
    return result

def triplets(size):
    """
    Returns a list of triplets with different indices
    in the range 0..size-1.
    """
    result = []
    for a in range(size):
        for b in range(a+1,size):
            for c in range(b+1,size):
               result.append((a, b, c))
    return result

def evaluate(powers, coeffs, point):
    """
    Evaluates the tropical polynomial at the point.
    """
    result = []
    for ind in range(len(coeffs)):
        (cff, pwr) = (coeffs[ind], powers[ind])
        val = cff + pwr[0]*point[0] + pwr[1]*point[1]
        result.append(val)
    return result

def minoccurstwice(values):
    """
    Returns true if the minimum occurs twice
    in the list of values.
    """
    minval = min(values)
    ind = values.index(minval)
    for k in range(len(values)):
        if(k != ind):
            dff = values[k] - minval 
            if(abs(dff) < TOL):
                return True
    return False

def triple(deg, cff, trp, verbose=True):
    """
    Computes the intersection of the equations
    defined by the indices in trp.
    Returns None if the intersection is not a minimum,
    otherwise returns a node, stored as a 2-tuple.
    The first element of the 2-tuple on return is the triplet trp,
    the second element is the tuple with the (x, y)-coordinates
    of the point where the minimum is attained three times.
    """
    result = None
    if verbose:
        print(deg[trp[0]], deg[trp[1]], deg[trp[2]])
    a = (deg[trp[1]][0] - deg[trp[0]][0], deg[trp[1]][1] - deg[trp[0]][1])
    b = (deg[trp[2]][0] - deg[trp[0]][0], deg[trp[2]][1] - deg[trp[0]][1])
    d = a[0]*b[1] - b[0]*a[1]
    if verbose:
        if(d == 0):
            print('trp = ', trp, 'det =', d)
    if(d != 0):
        r0 = cff[trp[0]] - cff[trp[1]]
        r1 = cff[trp[0]] - cff[trp[2]]
        x = (r0*b[1] - r1*a[1])/float(d)
        y = (a[0]*r1 - b[0]*r0)/float(d)
        if verbose:
            print('trp = ', trp, 'det =', d, x, y)
        val = cff[trp[0]] + deg[trp[0]][0]*x + deg[trp[0]][1]*y
        minval = min(evaluate(deg, cff, (x,y)))
        if verbose:
            print(cff[trp[0]] + deg[trp[0]][0]*x + deg[trp[0]][1]*y, \
                  cff[trp[1]] + deg[trp[1]][0]*x + deg[trp[1]][1]*y, \
                  cff[trp[2]] + deg[trp[2]][0]*x + deg[trp[2]][1]*y, minval)
        if(abs(val - minval) < TOL):
            if verbose:
                print('okay')
            result = (trp, (x, y))
    return result

def computenodes(powers, coeffs):
    """
    Computes all nodes of a tropical polynomial with
    exponents in powers and coefficients in coeffs.
    A node is where the minimum is attained at least 3 times.
    The list on return is a list of tuples, each tuples
    defines a node as a triplet of indices to the exponents
    and one 2-tuple with the coordinates of the node.
    """
    trp = triplets(len(coeffs))
    print('number of triplets :', len(trp))
    nodes = []
    for t in trp:
        node = triple(powers, coeffs, t, False)
        if(node != None):
            nodes.append(node)
    return nodes

def drawnodes(fig, axs, nodes, axsreset=True, axslim=None):
    """
    Draws dots to represent the nodes in the list.
    The list of four numbers axslim are used 
    in axs.set_xlim and axs.setylim.
    """
    xpoints = []
    ypoints = []
    for node in nodes:
        xpoints.append(node[1][0])
        ypoints.append(node[1][1])
    # print xpoints, ypoints
    if axsreset:
        if(axslim == None):
            axs.set_xlim(min(xpoints)-1.0, max(xpoints)+1.0)
            axs.set_ylim(min(ypoints)-1.0, max(ypoints)+1.0)
        else:
            axs.set_xlim(axslim[0], axslim[1])
            axs.set_ylim(axslim[2], axslim[3])
    dots, = axs.plot(xpoints, ypoints, 'ro')
    fig.canvas.draw()

def drawsegments(fig, axs, powers, coeffs, nodes, ind, nbr):
    """
    Draws the segments between the nodes in the list,
    in particular to nodes[ind] and all other nodes.
    The list nbr contains the number of connections made.
    """
    firstnode = nodes[ind]
    trpfirst = firstnode[0]
    xpt = firstnode[1][0]
    ypt = firstnode[1][1]
    for secondind in range(len(nodes)):
        if(secondind != ind):
            secondnode = nodes[secondind]
            trpsecond = secondnode[0]
            if((trpfirst[0] in trpsecond) or \
               (trpfirst[1] in trpsecond) or \
               (trpfirst[2] in trpsecond)):
                apt = secondnode[1][0]
                bpt = secondnode[1][1]
                xmd = (xpt + apt)/2.0
                ymd = (ypt + bpt)/2.0
                val = evaluate(powers, coeffs, (xmd, ymd))
                if minoccurstwice(val):
                    axs.plot([xpt, apt], [ypt, bpt], 'b-')
                    fig.canvas.draw()
                    nbr[0] = nbr[0] + 1

def drawrays(fig, axs, powers, coeffs, nodes, ind, it1, it2, nbr):
    """
    Draws the half rays out of the current node at nodes[ind],
    investigatings equations with indices it1 and it2 at the triplet.
    The list nbr contains the number of connections made.
    """
    firstnode = nodes[ind]
    trpfirst = firstnode[0]
    xpt = firstnode[1][0]
    ypt = firstnode[1][1]
    a = powers[trpfirst[it1]][0] - powers[trpfirst[it2]][0]
    b = powers[trpfirst[it1]][1] - powers[trpfirst[it2]][1]
    apt = xpt - OUTMAX*b
    bpt = ypt + OUTMAX*a
    val = evaluate(powers, coeffs, (apt, bpt))
    if minoccurstwice(val):
        dots, = axs.plot([xpt, apt], [ypt, bpt], 'r-')
        fig.canvas.draw()
        nbr[0] = nbr[0] + 1
    apt = xpt + OUTMAX*b
    bpt = ypt - OUTMAX*a
    val = evaluate(powers, coeffs, (apt, bpt))
    if minoccurstwice(val):
        dots, = axs.plot([xpt, apt], [ypt, bpt], 'r-')
        fig.canvas.draw()
        nbr[0] = nbr[0] + 1

def showcurve(fig, axs, nodes, powers, coeffs, axsset=True, axslim=None):
    """
    Plots the nodes in the list on the figure fig with axes axs,
    for the tropical curve with exponents in powers and 
    corresponding coeffficients in coeffs.
    """
    drawnodes(fig, axs, nodes, axsset, axslim)
    for firstind in range(len(nodes)):
        firstnode = nodes[firstind]
        xpt = firstnode[1][0]
        ypt = firstnode[1][1]
        trpfirst = firstnode[0]
        nbr = [0]
        drawsegments(fig, axs, powers, coeffs, nodes, firstind, nbr)
        drawrays(fig, axs, powers, coeffs, nodes, firstind, 0, 1, nbr)
        drawrays(fig, axs, powers, coeffs, nodes, firstind, 0, 2, nbr)
        drawrays(fig, axs, powers, coeffs, nodes, firstind, 1, 2, nbr)
        fig.canvas.draw()
        print('node', firstind, 'has', nbr[0], 'connections')

def drawbezout(fig, axs):
    """
    Illustrates the theorem of Bezout for two quadratic tropical curves.
    """
    from random import uniform as u
    deg = 2
    pwr = exponents(deg)
    pcff = []
    qcff = []
    # qcff = [u(-1, +1) for p in pwr]
    for expdeg in pwr:
        pcff.append(expdeg[0]**2 + expdeg[1]**2)
        qcff.append((3 - expdeg[0])**2 - (3 - expdeg[1])**2)
    pnodes = computenodes(pwr,pcff)
    showcurve(fig,axs,pnodes,pwr,pcff,True,[-10, 10,-15, 5])
    fig.canvas.draw()
    qnodes = computenodes(pwr,qcff)
    showcurve(fig,axs,qnodes,pwr,qcff,False)
    fig.canvas.draw()
    ans = input('hit enter to exit')

def drawstable(fig, axs):
    """
    Illustrates the stable intersection of a cubic with itself.
    """
    from random import uniform as u
    deg = 3
    pwr = exponents(deg)
    while True:
        pcff = [u(-1, +1) for p in pwr]
        pnodes = computenodes(pwr,pcff)
        showcurve(fig,axs,pnodes,pwr,pcff,True)
        fig.canvas.draw()
        # qnodes = nearbynodes(pnodes)
        qcff = []
        for pc in pcff:
            qcff.append(pc + u(+0.6,+0.8))
        qnodes = computenodes(pwr,qcff)
        showcurve(fig,axs,qnodes,pwr,qcff,False)
        fig.canvas.draw()
        ans = input('one other random case ? (y/n) ')
        if(ans != 'y'): break
        fig = plt.figure()
        axs = fig.add_subplot(111)

def drawcurve(fig, axs):
    """
    Prompt the user for a degree, generates random coefficients
    and then plots the tropical curve.
    """
    print('drawing a plane tropical curve ...')
    deg = input('give the degree : ')
    pwr = exponents(deg)
    print('the exponents :', pwr)
    prb = input('parabolic coefficients ? (y/n) ')
    from random import uniform as u
    while True:
        if(prb != 'y'):
            cff = [u(-1, +1) for p in pwr]
        else:
            cff = []
            for expdeg in pwr:
                cff.append(expdeg[0]**2 + expdeg[1]**2)
        print('the coefficients :', cff)
        nodes = computenodes(pwr, cff)
        print('the nodes :')
        for ind in range(len(nodes)):
            node = nodes[ind]
            print('node', ind, 'is', node)
        showcurve(fig, axs, nodes, pwr, cff)
        ans = input('draw for different coefficients ? (y/n) ')
        if(ans != 'y'): 
            break
        else:
            fig = plt.figure()
            axs = fig.add_subplot(111)

def main():
    """
    Prompts the user for a degree
    and generates random coefficients.
    """
    plt.ion()
    fig = plt.figure()
    axs = fig.add_subplot(111)
    ans = input("see an illustration of Bezout's theorem ? (y/n) ")
    if(ans == 'y'):
        drawbezout(fig, axs)
    else:
        ans = input("see an illustration of stable intersection ? (y/n) ")
        if(ans == 'y'):
            drawstable(fig, axs)
        else:
            drawcurve(fig, axs)

main()