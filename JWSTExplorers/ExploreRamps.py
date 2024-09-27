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
# from jwst.datamodels import dqflags
from jwst.datamodels.dqflags import group as grpdict
from jwst.datamodels.dqflags import pixel as pixeldict
# from stcal.dqflags import group as grpdict
from stcal.dqflags import interpret_bit_flags, dqflags_to_mnemonics

# sys.tracebacklimit = 0

class Explorer():
    def __init__(self,file):
        folder = os.path.dirname(file)+'/../'
        filename = os.path.basename(file)
        exten = '_'+filename.split('_')[-1].replace('.fits','')
        print(folder, filename, exten)
        print('\n')
        nfiles = 0
        if os.path.exists(folder+'stage0/'+filename.replace(exten,'_uncal')):
            hdu_uncal = fits.open(folder +'stage0/'+filename.replace(exten, '_uncal'))
            self.uncal_data = hdu_uncal['SCI'].data
            self.uncal = True
        else:
            self.uncal = False
            print('---- No _uncal file found ----')
        if os.path.exists(folder+'stage1/'+filename.replace(exten,'_fitopt')):
            self.hdu_fitopt = fits.open(folder+'stage1/'+filename.replace(exten,'_fitopt'))
            self.fitopt = True
        else:
            self.fitopt = False
            print('---- No _fitopt file found ----')
        if os.path.exists(folder+'stage1/'+filename.replace(exten,'_ramp')):
            hdu_ramp = fits.open(folder+'stage1/'+filename.replace(exten,'_ramp'))
            self.ramp_data = hdu_ramp['SCI'].data
            self.groupdq = hdu_ramp['GROUPDQ'].data
            self.pixeldq = hdu_ramp['PIXELDQ'].data
            self.img_data = self.ramp_data[0,-1,:,:]
            self.img_hdr0 = hdu_ramp[0].header
            self.img_hdr1 = hdu_ramp['SCI'].header
            self.ramps = True
        elif os.path.exists(folder+'stage1/'+filename.replace(exten,'_jump')):
            hdu_ramp = fits.open(folder+'stage1/'+filename.replace(exten,'_jump'))
            self.ramp_data = hdu_ramp['SCI'].data
            self.groupdq = hdu_ramp['GROUPDQ'].data
            self.pixeldq = hdu_ramp['PIXELDQ'].data
            self.img_data = self.ramp_data[0,-1,:,:]
            self.img_hdr0 = hdu_ramp[0].header
            self.img_hdr1 = hdu_ramp['SCI'].header
            self.ramps = True
        else:
            self.ramps = False
            print('---- No ramp file found ----')
        if os.path.exists(folder+'stage1/'+filename.replace(exten,'_rateints')):
            hdu_rate = fits.open(folder+'stage1/'+filename.replace(exten,'_rateints'))
            self.rate_data = hdu_rate['SCI'].data
            self.rate = True
        else:
            self.rate = False
            print('---- No _rateints file found ----')
        self.fig = plt.figure(figsize=(16,7))
        self.gs = gridspec.GridSpec(ncols=2,nrows=2,width_ratios=[1,1], height_ratios=[1,1])
        self.gs.update(wspace=0.15, hspace=0.3,left=0.1, bottom=0.1, right=0.9, top=0.9)
        self.ax1 = plt.subplot(self.gs[:,:1])
        self.ax2 = plt.subplot(self.gs[0,1])
        self.ax3 = plt.subplot(self.gs[1,1])
        self.divider = make_axes_locatable(self.ax1)
        self.cax = self.divider.append_axes("right", size="5%", pad=0.05)

    def start(self):
        print(
            '\n',
            '#' * 30,
        )
        print(
            'Note: Move your mouse over the image to plot the ramp'
        )
        print(
            'Note: Press the spacebar on the image to save the ramp of the corresponding pixel in the bottom-right subplot.'
        )
        print('#' * 30, '\n')

        # Initialize the left view
        scale = interval.ZScaleInterval(n_samples=600, contrast=0.25, max_reject=0.5, min_npixels=5, krej=2.5,max_iterations=5)
        (vmin, vmax) = scale.get_limits(self.img_data)
        normalization = matcol.Normalize(vmin=vmin, vmax=vmax)
        self.im = self.ax1.imshow(self.img_data, origin = 'lower', norm=normalization,cmap = cm.gray)
        cbar = plt.colorbar(self.im, cax=self.cax)
        cbar.set_label('Pixel unit (%s)'%self.img_hdr1['BUNIT'])
        self.ax1.set_xlabel('xpixel')
        self.ax1.set_ylabel('ypixel')
        self.ax1.set_title('%s\n%s'%(self.img_hdr0['TARGNAME'],self.img_hdr0['FILENAME']))
        bbox = dict(facecolor='white', alpha=0.7, edgecolor='grey', boxstyle='round')
        self.annot = self.ax1.annotate("", xy=(0,0), xytext=(10,10),textcoords="offset points", color='black', ha='center',bbox = bbox)
        self.patch1 = self.ax1.add_patch(plt.Rectangle((1e-5,1e-5), 1e-5,1e-5, ec = 'yellow',fc='None',linestyle='--', linewidth=1e-5))
        self.patch2 = self.ax1.add_patch(plt.Rectangle((1e-5,1e-5), 1e-5,1e-5, ec = 'yellow',fc='None',linestyle='--', linewidth=1e-5))
        self.patch3 = self.ax1.add_patch(plt.Rectangle((1e-5, 1e-5), 1e-5, 1e-5, ec='cyan', fc='None', linestyle='--',linewidth=1e-5))
        self.patch4 = self.ax1.add_patch(plt.Rectangle((1e-5, 1e-5), 1e-5, 1e-5, ec='cyan', fc='None', linestyle='--',linewidth=1e-5))
        # Initialize the commands
        self.mouseId = self.fig.canvas.mpl_connect("motion_notify_event", self.mouseMove)
        self.keyId = self.fig.canvas.mpl_connect("key_press_event", self.keyPress)
        plt.show()

    def mouseMove(self,event):
        self.annot.set_visible(False)
        self.patch1.set_visible(False)
        self.patch2.set_visible(False)
        self.ax2.clear()
        self.ax2.set_xlabel('Group number')
        self.ax2.set_ylabel('Count (DN)')
        # self.ax2.set_xlim([0.5,self.ramp_data.shape[1]*5+0.5])
        self.ax2.set_xlim([0.5,self.ramp_data.shape[1]*self.ramp_data.shape[0]+0.5])
        try:
            origin = np.array((round(event.xdata), round(event.ydata)))
            if np.all(origin>=0) and np.all(origin<=max(self.ramp_data.shape[-2],self.ramp_data.shape[-1])):
                for j in range(self.ramp_data.shape[0]):
                    if self.uncal and self.ramps:
                        self.ax2.plot([(i+self.uncal_data.shape[1]*j)+1 for i in range(self.uncal_data.shape[1])], [self.uncal_data[j][i][origin[1]][origin[0]] for i in range(len(self.uncal_data[j]))],marker='s',markerfacecolor='lightgrey',markeredgecolor='lightgrey', markersize=4,linestyle='',alpha=0.8)
                    self.ax2.plot([(i+self.ramp_data.shape[1]*j)+1 for i in range(self.ramp_data.shape[1])], [self.ramp_data[j][i][origin[1]][origin[0]] for i in range(len(self.ramp_data[j]))],marker='o',markerfacecolor='#0098FF',markeredgecolor='#0098FF', markersize=4,linestyle='')
                    if self.fitopt:
                        slope = self.hdu_fitopt['SLOPE'].data[j,:,origin[1],origin[0]]
                        zp = self.hdu_fitopt['YINT'].data[j,:,origin[1],origin[0]]
                        var_p = self.hdu_fitopt['VAR_POISSON'].data[j,:,origin[1],origin[0]]
                        var_r = self.hdu_fitopt['VAR_RNOISE'].data[j,:,origin[1],origin[0]]
                        var_comb = var_r + var_p
                        var_comb[var_comb==0] = np.nan
                        if np.nansum(var_comb)==0:
                            pass
                        else:
                            slope_int = np.nansum(np.nan_to_num(slope/var_comb, posinf=0, nan=np.nan)) / np.nansum(np.nan_to_num(1/var_comb, posinf=0, nan=np.nan))
                            zp_int = np.nansum(np.nan_to_num(zp/var_comb, posinf=0, nan=np.nan)) / np.nansum(np.nan_to_num(1/var_comb, posinf=0, nan=np.nan)) # np.sum(np.nan_to_num(1/var_comb,posinf=0))
                            t1,t2 = -self.img_hdr0['TGROUP'], (self.ramp_data.shape[1]-3)*self.img_hdr0['TGROUP']
                            self.ax2.plot([1+self.ramp_data.shape[1]*j,self.ramp_data.shape[1]+self.ramp_data.shape[1]*j-1], [zp_int+t1*slope_int, zp_int+t2*slope_int],linestyle='-',color='black',label='Slope %i = %.2f DN/s'%(j+1,slope_int))
                    elif self.rate:
                        t1, t2 = -self.img_hdr0['TGROUP'], (self.ramp_data.shape[1]-3)*self.img_hdr0['TGROUP']
                        zp_int = self.ramp_data[0,1,origin[1],origin[0]]
                        slope_int = self.rate_data[0,origin[1],origin[0]]
                        self.ax2.plot([1+self.ramp_data.shape[1]*j,self.ramp_data.shape[1]+self.ramp_data.shape[1]*j-1], [zp_int+t1*slope_int, zp_int+t2*slope_int],linestyle='-',color='black',label='Slope %i = %.2f DN/s'%(j+1,slope_int))
                self.ramp(self.ax2, origin[0], origin[1])
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
                    pxdqvalue = self.pixeldq[origin[1]][origin[0]]
                    pixflag = dqflags_to_mnemonics(pxdqvalue, pixeldict)
                    if len(pixflag)==0:
                        self.annot.set_text("OK")
                    else:
                        self.annot.set_text("%s: %s"%(pxdqvalue, list(pixflag)))
                    self.annot.set_visible(True)
                square = plt.Rectangle((origin[0] - 0.5, origin[1] - 0.5), 1, 1, ec='yellow', fc='None')
                circle = plt.Circle((origin[0], origin[1]), 0.02 * self.ramp_data.shape[-1]+self.ramp_data.shape[-2]/self.ramp_data.shape[-1], ec='yellow', fc='None',
                                     linestyle='--', linewidth=0.6)
                self.patch1 = self.ax1.add_patch(square)
                self.patch2 = self.ax1.add_patch(circle)
                self.patch1.set_visible(True)
                self.patch2.set_visible(True)
                self.ax2.set_title('Pixel %.2fx%.2f (%ix%i)'%(event.xdata,event.ydata,origin[0],origin[1]))
                self.ax2.set_title('Pixel %ix%i'%(origin[0],origin[1]))
        except:
            pass
        plt.draw()

    def keyPress(self, event):
        sys.stdout.flush()
        if event.key == ' ':
            self.ax3.clear()
            self.patch3.remove()
            self.patch4.remove()
            origin = (round(event.xdata), round(event.ydata))
            print('press spacebar: x = %i, y = %i'%(origin[0], origin[1]))
            square = plt.Rectangle((origin[0]-0.5,origin[1]-0.5), 1,1, ec = 'orangered',fc='None')
            circle = plt.Circle((origin[0],origin[1]), 0.02*self.ramp_data.shape[-1]+self.ramp_data.shape[-2]/self.ramp_data.shape[-1], ec = 'orangered',fc='None',linestyle='--', linewidth=0.6)
            self.patch3 = self.ax1.add_patch(square)
            self.patch4 = self.ax1.add_patch(circle)

            for j in range(self.ramp_data.shape[0]):
                if self.uncal and self.ramps:
                    self.ax3.plot([(i+self.uncal_data.shape[1]*j)+1 for i in range(self.uncal_data.shape[1])], [self.uncal_data[j][i][origin[1]][origin[0]] for i in range(len(self.uncal_data[j]))],marker='s',markerfacecolor='lightgrey',markeredgecolor='lightgrey', markersize=4,linestyle='',alpha=0.8)
                self.ax3.plot([(i+self.ramp_data.shape[1]*j)+1 for i in range(self.ramp_data.shape[1])], [self.ramp_data[j][i][origin[1]][origin[0]] for i in range(len(self.ramp_data[j]))],marker='o',markerfacecolor='#0098FF',markeredgecolor='#0098FF', markersize=4,linestyle='')
                if self.fitopt:
                    slope = self.hdu_fitopt['SLOPE'].data[j,:,origin[1],origin[0]]
                    zp = self.hdu_fitopt['YINT'].data[j,:,origin[1],origin[0]]
                    var_p = self.hdu_fitopt['VAR_POISSON'].data[j,:,origin[1],origin[0]]
                    var_r = self.hdu_fitopt['VAR_RNOISE'].data[j,:,origin[1],origin[0]]
                    var_comb = var_r + var_p
                    var_comb[var_comb==0] = np.nan
                    slope_int = np.sum(np.nan_to_num(slope/var_comb))/np.sum(np.nan_to_num(1/var_comb,posinf=0))
                    zp_int = np.sum(np.nan_to_num(zp/var_comb))/np.sum(np.nan_to_num(1/var_comb,posinf=0))
                    # zp = self.hdu2['YINT'].data[j][origin[1]][origin[0]]
                    # slope = self.hdu2['SLOPE'].data[j][origin[1]][origin[0]]
                    t1,t2 = -self.img_hdr0['TGROUP'], (self.ramp_data.shape[1]-3)*self.img_hdr0['TGROUP']
                    self.ax3.plot([1+self.ramp_data.shape[1]*j,self.ramp_data.shape[1]+self.ramp_data.shape[1]*j-1], [zp_int+t1*slope_int, zp_int+t2*slope_int],linestyle='-',color='black',label='Slope %i = %.2f DN/s'%(j+1,slope_int))
                elif self.rate:
                    t1, t2 = -self.img_hdr0['TGROUP'], (self.ramp_data.shape[1]-3)*self.img_hdr0['TGROUP']
                    zp_int = self.ramp_data[0,1,origin[1],origin[0]]
                    slope_int = self.rate_data[0,origin[1],origin[0]]
                    self.ax3.plot([1+self.ramp_data.shape[1]*j,self.ramp_data.shape[1]+self.ramp_data.shape[1]*j-1], [zp_int+t1*slope_int, zp_int+t2*slope_int],linestyle='-',color='black',label='Slope %i = %.2f DN/s'%(j+1,slope_int))
            self.ramp(self.ax3, origin[0],origin[1])
            self.ax3.set_ylabel('Count (DN)')
            self.ax3.set_xlabel('Group number')
            # self.ax3.set_xlim([0.5,self.ramp_data.shape[1]*5+0.5])
            self.ax3.set_xlim([0.5,self.ramp_data.shape[1]*self.ramp_data.shape[0]+0.5])
            self.ax3.set_title('Pixel %ix%i'%(origin[0],origin[1]))
            plt.draw()

    def ramp(self,ax,xpixel,ypixel):
        flag_d, flag_j, flag_s,flag_e = -1, -1, -1, -1
        for j in range(self.groupdq.shape[0]):
            for i in range(len(self.groupdq[j])):  
                flag = dqflags_to_mnemonics(self.groupdq[j][i][ypixel][xpixel], grpdict)
                flag = list(flag)
                length = 0
                if 'JUMP_DET' in flag:
                    length+=1
                    flag_j = i
                    ax.plot(i+1+self.ramp_data.shape[1]*j, self.ramp_data[j][i][ypixel][xpixel],marker='o',markerfacecolor="None", markeredgecolor='red', label=None,markersize=4,linestyle='')
                if 'SATURATED' in flag:
                    length+=1
                    flag_s = i
                    ax.plot(i+1+self.ramp_data.shape[1]*j, self.ramp_data[j][i][ypixel][xpixel],marker='o',markerfacecolor="yellow", markeredgecolor='yellow', label=None,markersize=4,linestyle='')
                if 'DO_NOT_USE' in flag:
                    length+=1
                    flag_d = i
                    ax.plot(i+1+self.ramp_data.shape[1]*j, self.ramp_data[j][i][ypixel][xpixel],marker='x',markerfacecolor="None", markeredgecolor='black', label=None,markersize=4,linestyle='')
                if length!=len(flag):
                    flag_e = i
                    ax.plot(i+1+self.ramp_data.shape[1]*j, self.ramp_data[j][i][ypixel][xpixel],marker='+',markerfacecolor="None", markeredgecolor='green', label=None,markersize=4,linestyle='')
        if flag_j!=-1:
            ax.plot([], [],marker='o',markerfacecolor="None",markeredgecolor='red', markersize=4, label='JUMP_DET',linestyle='')
        if flag_s!=-1:
            ax.plot([], [],marker='o',markerfacecolor="yellow",markeredgecolor='yellow', markersize=4, label='SATURATED',linestyle='')
        if flag_d!=-1:
            ax.plot([], [],marker='x',markerfacecolor="None",markeredgecolor='black', markersize=4, label='DO_NOT_USE',linestyle='')
        if flag_e!=-1:
            ax.plot([], [],marker='+',markerfacecolor="None", markeredgecolor='green', label='OTHER FLAG',markersize=4,linestyle='')
        if self.uncal:
            ax.plot([], [], marker='s', markerfacecolor='lightgrey', markeredgecolor='lightgrey', markersize=4, label='Uncalibrated', linestyle='')
        ax.plot([],[], marker='o', markerfacecolor='lightgrey', markeredgecolor='lightgrey', markersize=4, label='Calibrated', linestyle='')
        ax.legend(loc='upper left', fontsize = 8, ncol=3)

if __name__ == "__main__":
    explorer = Explorer(str(sys.argv[1]))
    explorer.start()