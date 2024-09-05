import pylab as pl
from milkywayplots.milkywayplot import make_mw_plot, get_image

mw_img_name = get_image()

fig = pl.figure(figsize=(8,8), num=5)
fig.clf()
ax, ax_pixgrid, gcp, hcp, tr_gal, tr_helio = make_mw_plot(mw_img_name=mw_img_name)
hcp.plot(0,0,'o',color='gold',markersize=10,markeredgecolor='gold',markerfacecolor='none',zorder=50,alpha=1,markeredgewidth=2)
hcp.plot(0,0,'.',color='gold',markersize=3,markeredgecolor='gold',markerfacecolor='gold',zorder=50,alpha=1,markeredgewidth=2)

# GC
gcp.plot(0, 0, 'o', color='red', markersize=40, markerfacecolor='none', zorder=50, alpha=1, markeredgewidth=2, markeredgecolor='r')
gcp.text(270, 1.5, 'Galactic Center', ha='center', color='k')

# wd 2
hcp.plot(284, 4.03, 'o', color='red', markersize=10, markerfacecolor='none', zorder=50, alpha=1, markeredgewidth=2, markeredgecolor='r')
hcp.text(278, 4.03, 'Wd 2', color='white', ha='center', va='bottom')
# wd 1
hcp.plot(339.54, 3.8, 'o', color='red', markersize=10, markerfacecolor='none', zorder=50, alpha=1, markeredgewidth=2, markeredgecolor='r')
hcp.text(334, 3.6, 'Wd 1', color='white', ha='center', va='bottom')

hcp.plot(291.6240, 6.9, 'o', color='red', markersize=10, markerfacecolor='none', zorder=50, alpha=1, markeredgewidth=2, markeredgecolor='r')
hcp.text(287, 6.9, 'NGC 3603', color='white', ha='center', va='bottom')

hcp.plot(49.5, 5.1, 'o', color='green', markersize=10, markerfacecolor='none', zorder=50, alpha=1, markeredgewidth=2, markeredgecolor='green')
hcp.text(54, 5.1, 'W51', color='white', ha='center', va='bottom')

hcp.plot(30, 6.0, 'o', color='green', markersize=10, markerfacecolor='none', zorder=50, alpha=1, markeredgewidth=2, markeredgecolor='green')
hcp.text(34, 6.0, 'W43', color='white', ha='center', va='bottom')

ax.axis['bottom'].major_ticklabels.set_visible(False)
ax.axis['left'].major_ticklabels.set_visible(False)


pl.draw()
pl.savefig("galaxyoverview_clusters_and_protoclusters.png", bbox_inches='tight')
ax.axis([-10, 10, -10, 10])
pl.savefig("galaxyoverview_clusters_and_protoclusters_zoom.png", bbox_inches='tight')

for elt in hcp.texts + gcp.texts + gcp.lines + hcp.lines:
    elt.set_visible(False)

hcp.plot(0,0,'o',color='gold',markersize=10,markeredgecolor='gold',markerfacecolor='none',zorder=50,alpha=1,markeredgewidth=2)
hcp.plot(0,0,'.',color='gold',markersize=3,markeredgecolor='gold',markerfacecolor='gold',zorder=50,alpha=1,markeredgewidth=2)

# ww, hh = 1, 1
# rect = pl.matplotlib.patches.Rectangle([0-ww/2, 8.5-hh/2], ww, hh)
# rect.set_edgecolor('red')
# rect.set_facecolor('none')
# ax.add_artist(rect)

circ = pl.matplotlib.patches.Circle([0, 8.5], 1)
circ.set_edgecolor('red')
circ.set_facecolor((1, 0, 0, 0.25))
ax.add_artist(circ)
pl.savefig("galaxyoverview_local_kiloparsec.png", bbox_inches='tight', edgecolor='none', facecolor='none')
