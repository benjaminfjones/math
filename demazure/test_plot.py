# plot styles 
arrow_style={'rgbcolor':(0,0,1),'width':1,'linestyle':'dotted'}
point_style={'rgbcolor':(1,0,0)}
grid_style = {'alpha':0.5, 'rgbcolor':(0,1,0), 'linestyle':'--'}
mult_scale = lambda m: 8 * (m + 1)**2

# roots and weights
ct = ('G',2)
wcr = WeylCharacterRing(ct)
wr = WeightRing(wcr)
space = wcr.space()
rts = space.roots()
vroots = [ vector(r) for r in rts ]
coroots = [ (2/norm(r)**2)*r for r in vroots ]
scoroots = [ (1/norm(r)**2)*r for r in coroots ]
lam = space.fundamental_weights()

# module
M1 = wcr(lam[1])
M2 = wcr(lam[2])
ch = WeightRingElement(wr, M1._mdict) + WeightRingElement(wr, M2._mdict)

root_coords = [vector((r[1], r[0])) for r in roots]
scoroot_coords = [vector((r[1], r[0])) for r in scoroots]
wt_coords = [(vector((wt[1], wt[0])), m) for (wt, m) in ch.mlist()]
lam_coords = [vector((l[1], l[0])) for l in lam]

# root and weight plots
root_plot = sum(plot(v,**arrow_style) for v in root_coords)
coroot_plot = sum(plot(v,rgbcolor=(1,0,1)) for v in scoroot_coords)
wt_plot = sum(point(wt, size=mult_scale(m), **point_style) 
                       for (wt, m) in wt_coords)
roots_and_weights = root_plot + wt_plot
mmdata = roots_and_weights.get_minmax_data()
xm = max(abs(mmdata['xmax']), abs(mmdata['xmin']))
ym = max(abs(mmdata['ymax']), abs(mmdata['ymin']))
d = ceil(sqrt(xm*xm + ym*ym))

grid_lines = []
for v in scoroot_coords:
    # non-vertical case
    if v[1] != 0:
        s = QQ(-v[0])/QQ(v[1]) # slope
        p1 = vector((d, s*d))
    # vertical 
    else:
        p1 = vector((0,d))
    for k in range(-d,d+1):
        # plot the k-shifted hyperplane corresponding to l
        grid_lines.append(line([p1 + k*v, -p1 + k*v], **grid_style))
#   grid_lines.append(plot(lv, rgbcolor=(1,0,1)))

#show(roots_and_weights, aspect_ratio=1, axes=False, **mmdata)
show(sum(grid_lines) + roots_and_weights, aspect_ratio=1, axes=False, **mmdata)