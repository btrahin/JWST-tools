import sys, os, glob
import numpy as np
import matplotlib
# matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as matcol
from matplotlib import cm, gridspec
import mpl_interactions.ipyplot as iplt #TO DO
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from itertools import combinations
from astropy.visualization import interval
from astropy.io import fits
from jwst.datamodels.dqflags import pixel as pixeldict
from stcal.dqflags import interpret_bit_flags, dqflags_to_mnemonics

# sys.tracebacklimit = 0


class Explorer():

    def __init__(self,file):
        print('Working on %s'%file)
        print('\n')
        with fits.open(file) as hdu:
            self.head0 = hdu[0].header
            self.head1 = hdu[1].header
            self.data = hdu[1].data
            self.flag = hdu['DQ'].data

        self.fig = plt.figure(figsize=(14,7))
        self.gs = gridspec.GridSpec(ncols=2,nrows=1,width_ratios=[1,1], height_ratios=[1])
        self.gs.update(wspace=0.3, hspace=0.3,left=0.06, bottom=0.08, right=0.99, top=0.95)
        self.ax1 = plt.subplot(self.gs[0,0])
        self.ax2 = plt.subplot(self.gs[0,1], sharex=self.ax1, sharey=self.ax1)
        self.divider = make_axes_locatable(self.ax1)
        self.cax = self.divider.append_axes("right", size="5%", pad=0.05)

    def start(self):
        print(
            '\n',
            '#' * 30,
        )
        print(
            'Note: Move your mouse on the image to see the flag'
        )
        print(
            'Note: Click/drag your mouse on the image to flag a pixel as DO_NOT_USE'
        )
        print(
            'Note: Press c and drag on the image to create a circle of DO_NOT_USE pixels'
        )
        print(
            'Note: Press r and click/drag on the image to create a rectangle of DO_NOT_USE pixels'
        )
        print(
            'Note: Press s to save the fits file'
        )
        print('#' * 30, '\n')

        # Initialize the views
        scale = interval.ZScaleInterval(n_samples=600, contrast=0.25, max_reject=0.5, min_npixels=5, krej=2.5,max_iterations=5)
        self.image = self.data[:,:]
        (vmin, vmax) = scale.get_limits(self.image)
        normalization = matcol.Normalize(vmin=vmin, vmax=vmax)
        self.im = self.ax1.imshow(self.image,origin = 'lower', norm=normalization,cmap = cm.gray)
        cbar = plt.colorbar(self.im, cax=self.cax)
        cbar.set_label('Pixel unit (%s)'%self.head1['BUNIT'])
        self.ax1.set_xlabel('xpixel')
        self.ax1.set_ylabel('ypixel')
        self.ax1.set_title('%s\n%s'%(self.head0['TARGNAME'],self.head0['FILENAME']))
        self.image2 = self.flag[:,:]
        (vmin, vmax) = scale.get_limits(self.image2)
        normalization = matcol.Normalize(vmin=vmin, vmax=vmax)
        self.im2 = self.ax2.imshow(self.image2,origin = 'lower', norm=normalization,cmap = cm.gray)
        self.ax2.set_xlabel('xpixel')
        self.ax2.set_ylabel('ypixel')
        self.ax2.set_title('DQ image')
        # Initialize the commands
        bbox = dict(facecolor='white', alpha=0.7, edgecolor='grey', boxstyle='round')
        self.annot = self.ax1.annotate("", xy=(0,0), xytext=(10,10),textcoords="offset points", color='black', ha='center',bbox = bbox)
        self.patch1 = self.ax1.add_patch(plt.Rectangle((0.000001, 0.000001), 0.000001, 0.000001, ec='black', fc='None', linestyle='--',linewidth=0.000001))
        self.patch2 = self.ax1.add_patch(plt.Rectangle((0.000001, 0.000001), 0.000001, 0.000001, ec='grey', fc='white', alpha=0.4, linestyle='--',linewidth=0.000001))
        # self.annot2 = self.ax2.annotate("", xy=(0, 0), xytext=(10, 10), textcoords="offset points", color='black',
        #                                ha='center', bbox=bbox)
#         self.patch3 = self.ax2.add_patch(plt.Rectangle((0.000001, 0.000001), 0.000001, 0.000001, ec='yellow', fc='None', linestyle='--',linewidth=0.000001))
#         self.patch4 = self.ax2.add_patch(plt.Rectangle((0.000001, 0.000001), 0.000001, 0.000001, ec='yellow', fc='None', linestyle='--',linewidth=0.000001))
        self.mouseId = self.fig.canvas.mpl_connect("motion_notify_event", self.mouseMove)
        # self.keyPressId = self.fig.canvas.mpl_connect("key_press_event", lambda event: self.keyPress(event))
        # self.keyReleaseId = self.fig.canvas.mpl_connect("key_press_event", lambda event: self.keyRelease(event))
        plt.show()

    def mouseMove(self,event):
        self.annot.set_visible(False)
        self.patch1.set_visible(False)
        self.patch2.set_visible(False)
        # self.annot2.set_visible(False)
        # self.patch3.set_visible(False)
        # self.patch4.set_visible(False)
        try:
            origin = (round(event.xdata), round(event.ydata))
            if event.inaxes == self.ax1:
                self.annot.xy = origin
#                 self.annot2.xy = origin
                if origin[0]<=300:
                    self.annot.set_ha('left')
#                     self.annot2.set_ha('left')
                if 300<origin[0]<700:
                    self.annot.set_ha('center')
#                     self.annot2.set_ha('center')
                if 700<=origin[0]:
                    self.annot.set_ha('right')
#                     self.annot2.set_ha('right')
                dqvalue = self.flag[origin[1]][origin[0]]
                pixflag = dqflags_to_mnemonics(dqvalue, pixeldict)
                if len(pixflag)==0:
                    self.annot.set_text("NORMAL")
#                     self.annot2.set_text("NORMAL")
                else:
                    self.annot.set_text("%s: %s"%(dqvalue, list(pixflag)))
#                     self.annot2.set_text("%s: %s" % (dqvalue, list(pixflag)))
                self.annot.set_visible(True)
#                 self.annot2.set_visible(True)
                square = plt.Rectangle((origin[0] - 0.5, origin[1] - 0.5), 1, 1, ec='black', fc='None')
                circle = plt.Circle((origin[0], origin[1]), 0.02 * self.data.shape[-1], ec='grey', fc='white', alpha=0.4,
                                    linestyle='--', linewidth=0.6)
                self.patch1 = self.ax1.add_patch(square)
                self.patch2 = self.ax1.add_patch(circle)
                self.patch1.set_visible(True)
                self.patch2.set_visible(True)
                # self.patch3 = self.ax2.add_patch(square)
                # self.patch4 = self.ax2.add_patch(circle)
                # self.patch3.set_visible(True)
                # self.patch4.set_visible(True)
                plt.draw()
            else:
                self.annot.set_visible(False)
                self.patch1.set_visible(False)
                self.patch2.set_visible(False)
                # self.annot2.set_visible(False)
                # self.patch3.set_visible(False)
                # self.patch4.set_visible(False)
        except:
            pass

    def keyPress(self, event):
        if event.key == 'a':
            self.drawkey='pixel'
            self.mousePressId = self.fig.canvas.mpl_connect("button_press_event", lambda event: self.mousePress(event))
            self.mouseReleaseId = self.fig.canvas.mpl_connect("button_release_event", lambda event: self.mouseRelease(event))
        # elif event.key == "r":
        #     self.drawkey='rectangle'
        #     self.mousePressId = self.fig.canvas.mpl_connect("button_press_event", lambda event: self.mousePress(event))
        #     self.mouseReleaseId = self.fig.canvas.mpl_connect("button_release_event", lambda event: self.mouseRelease(event))
        # elif event.key == "c":
        #     self.drawkey='circle'
        #     self.mousePressId = self.fig.canvas.mpl_connect("button_press_event", lambda event: self.mousePress(event))
        #     self.mouseReleaseId = self.fig.canvas.mpl_connect("button_release_event", lambda event: self.mouseRelease(event))

    def keyRelease(self, event):
        pass

    def mousePress(self,event):
        self.mouseMouveId = self.fig.canvas.mpl_connect("motion_notify_event", lambda event: self.donotuse(event))

    def mouseRelease(self, event):
        self.fig.canvas.mpl_disconnect(self.mouseMouveId)

    def donotuse(self,event):
        if self.drawkey=='pixel':
            square = plt.Rectangle(
                (round(event.xdata) - 0.5, round(event.ydata) - 0.5), 1, 1, fc='yellow',alpha=0.7)
            self.ax1.add_patch(square)
            self.flag[round(event.ydata),round(event.xdata)]=1

    # def draw_click(self, event):
    #     # size = square (4 * duration of the time button 
    #     # is keep pressed )
    #     size = 4 * (self.end_time - self.start_time) ** 2
          
    #     # create a point of size=0.002 where mouse button 
    #     # clicked on the plot
    #     c1 = plt.Circle([event.xdata, event.ydata], 0.002,)
          
    #     # create a circle of radius 0.02*size
    #     c2 = plt.Circle([event.xdata, event.ydata], 0.02 * size, alpha=0.2)
    #     event.canvas.figure.gca().add_artist(c1)
    #     event.canvas.figure.gca().add_artist(c2)
    #     event.canvas.figure.show()

if __name__ == "__main__":
    explorer = Explorer(str(sys.argv[1]))
    explorer.start()