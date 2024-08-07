# python explore_cube.py to run
from os import EX_CANTCREAT
import sys
import glob
import numpy as np
import random
import matplotlib
# matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from matplotlib import cm, gridspec
import mpl_interactions.ipyplot as iplt
from astropy.io import fits
import warnings
import matplotlib.colors as matcol
from astropy.visualization import interval

sys.tracebacklimit = 0
# warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)


class Explorer():
    def __init__(self, filename):
        self.file = filename
        self.hdu = fits.open(self.file)
        self.fig = plt.figure(figsize=(14, 7))
        self.gs = gridspec.GridSpec(ncols=2,
                                    nrows=2,
                                    width_ratios=[1, 1],
                                    height_ratios=[1, 1])
        self.gs.update(hspace=0.4,
                       left=0.07,
                       bottom=0.08,
                       right=0.99,
                       top=0.95)
        self.ax1 = plt.subplot(self.gs[0, 0])
        self.ax2 = plt.subplot(self.gs[0, 1:])
        self.ax3 = plt.subplot(self.gs[1, 0])
        self.ax4 = plt.subplot(self.gs[1, 1])
        try: 
            self.data = self.hdu[1].data
            self.head = self.hdu[1].header
            self.pixscale = self.head['CDELT2']*3600 # deg in arcsec
            self.unit = self.head['BUNIT'] # MJy/sr
        except:
            self.data = self.hdu[0].data
            self.head = self.hdu[0].header
            self.pixscale = self.head['CDELT2'] # in arcsec
            self.unit = self.head['UNITS'] # muJy/arcsec2
        self.ymin, self.ymax = 0, np.nanmax(self.data[-1])#np.nanmax(self.data)
        self.nframes, self.nx, self.ny = self.data.shape
        try:
            self.wave = [
            self.head['CRVAL3'] + i * self.head['CDELT3']
            for i in range(self.head['NAXIS3'])
            ]
        except:
            self.wave = [i[0] for i in self.hdu['WCS-TABLE'].data[0][0]]
        self.ax2.set_xlabel(r'Wavelength ($\mu$m)')
        self.ax2.set_ylabel('Surface brightness (%s)'%self.unit)
        self.ax3.set_title(r'Average spectrum', fontsize=10)
        self.ax3.set_xlabel(r'Wavelength ($\mu$m)')
        self.ax3.set_ylabel('Integrated flux (Jy)')
        self.ax4.set_title(r'Saved spectra', fontsize=10)
        self.ax4.set_xlabel(r'Wavelength ($\mu$m)')
        self.ax4.set_ylabel('Surface brightness (%s)'%self.unit)

    def start(self):
        '''Initialize the cube and different panels
        '''
        print(
            '\n',
            '#' * 30,
        )
        print(
            'Note: Move your mouse on the cube to see the spectrum at each spaxel.'
        )
        print(
            'Note: You can see the different slices of the cube using the scroll bar on the second window'
        )
        print(
            'Note: Press x on the cube to save a spectrum on the subplot below.'
        )
        print('Note: Press c on the cube to erase all saved spectra.')
        print('Note: Press p on the cube to save spectrum extracted in a region of the size of the PSF in ch3 (1" diameter)')
        print('Note: Click and drag your mouse on the cube to save spectrum extracted in the region of the corresponding area')
        print(
            'Note: Press < on a saved spectrum to execute a function.'
        )
        print('#' * 30, '\n')
        # Initialize the left view
        # scale = interval.ZScaleInterval(nsamples=800, contrast=0.3, max_reject=0.5, min_npixels=5, krej=2.5,max_iterations=5)
        # (vmin, vmax) = scale.get_limits(self.data[:,:,:])
        # normalization = matcol.Normalize(vmin=vmin, vmax=vmax)
        self.im = iplt.imshow(lambda lambda_um: self.data[
            self.wave.index(lambda_um), :, :],
                              lambda_um=self.wave,
                              origin='lower',
                              # norm = normalization,
                              cmap=cm.viridis,
                              ax=self.ax1,
                              label="Wavelength (microns)")
        plt.colorbar()
        self.ax1.set_xlabel('xpixel')
        self.ax1.set_ylabel('ypixel')
        # Initialize the commands
        self.mouse = self.fig.canvas.mpl_connect("motion_notify_event",
                                                 self.mouseMove)
        self.key = self.fig.canvas.mpl_connect("key_press_event",
                                               self.keyPress)
        self.mousePressId = self.fig.canvas.mpl_connect("button_press_event", lambda event: self.mousePress(event))
        self.mouseReleaseId = self.fig.canvas.mpl_connect("button_release_event", lambda event: self.mouseRelease(event))
        
        # Initialize the rectangle
        self.rect = plt.Rectangle((0, 0), 0, 0, linewidth=1, edgecolor='grey', facecolor='grey', alpha=0.5)
        self.ax1.add_patch(self.rect)

        plt.show()

    def mousePress(self,event):
        '''Plot the average spectrum in the region defined in the cube
        '''
        origin = (round(event.xdata), round(event.ydata))
        
        # Update rectangle position
        self.rect.set_y(origin[1])
        self.rect.set_x(origin[0])

        # Activate the rectangle update
        self.mouseMouveId = self.fig.canvas.mpl_connect("motion_notify_event", lambda event: self.AverageSpectrum(event, origin))

    def mouseRelease(self, event):
        self.fig.canvas.mpl_disconnect(self.mouseMouveId)

    def AverageSpectrum(self, event, origin):
        '''Plot average spectrum (to be improved) in the lower left panel in the region defined by clicking and moving the mouse in the cube
        '''

        self.ax3.clear()
        self.ax3.set_xlabel(r'Wavelength ($\mu$m)')
        self.ax3.set_ylabel('Integrated flux (Jy)')
        # self.ax3.set_ylim(self.ymin, self.ymax)
        # self.ax3.set_xticks(np.arange(self.wave[0], self.wave[-1], 0.2))
        
        try:
            limits = np.array((round(event.xdata), round(event.ydata)))
            if np.all(limits > 0) and np.all(limits < self.nx):
                delta_x, delta_y = round(origin[0] - event.xdata), round(origin[1] - event.ydata)
                self.rect.set_width(-delta_x)
                self.rect.set_height(-delta_y)
                x_0 = min(origin[0], round(event.xdata))
                x_f = max(origin[0], round(event.xdata))
                y_0 = min(origin[1], round(event.ydata))
                y_f = max(origin[1], round(event.ydata))
                flux = self.data[:,y_0:y_f+1,x_0:x_f+1]*self.pixscale**2
                if 'sr' in self.unit: 
                    flux /= 4.255e10 # flux in MJy
                    flux /= 1e6 # flux in Jy
                else: 
                    flux *= 1e-6 # flux in Jy
                mean_flux = np.zeros(flux.shape[0])
                for i in range(flux.shape[0]):
                    mean_flux[i] = np.mean(flux[i])
                self.ax3.plot(self.wave,
                    mean_flux,
                    color='black',
                    linewidth=0.8)
                plt.draw()
        except:
            pass

    def mouseMove(self, event):
        '''Plot spectrum in the upper right panel by moving the mouse in the different spaxels of the cube
        '''
        self.ax2.clear()
        self.ax2.set_ylim(self.ymin, self.ymax)
        # self.ax2.set_xticks(np.arange(self.wave[0], self.wave[-1], 0.2))
        try:
            origin = np.array((round(event.xdata), round(event.ydata)))
            if np.all(origin > 0) and np.all(origin < self.nx):
                self.ax2.plot(self.wave, [
                    self.data[i][origin[1]][origin[0]]
                    for i in range(self.data.shape[0])
                ],
                              color='black',
                              linewidth=0.8)
                self.ax2.set_title('Pixel %ix%i' %
                                   (round(event.xdata), round(event.ydata)))
                plt.draw()
        except:
            pass

    def keyPress(self, event):
        '''You can assign here different keys to run functions inside the upper left panel
        '''
        sys.stdout.flush()
        if event.key == 'c':
            print('Spectra have been erased', event.key)
            self.ax4.clear()
            self.ax4.set_title(r'Saved spectra', fontsize=10)
            self.ax4.set_xlabel(r'Wavelength ($\mu$m)')
            self.ax4.set_ylabel('Surface brightness (%s)'%self.unit)
            plt.draw()
        elif event.key == 'x':
            self.origin = (round(event.xdata), round(event.ydata))
            print('Spectrum has been saved', event.key)
            print('x = %i, y = %i' % (self.origin[0], self.origin[1]))
            color = [random.random(), random.random(), random.random()]
            square = plt.Rectangle(
                (self.origin[0] - 0.5, self.origin[1] - 0.5), 1, 1, fc=color)
            self.ax1.add_patch(square)
            self.ax1.text(self.origin[0]+1, self.origin[1]+1, s="(%i, %i)"%(self.origin[0], self.origin[1]), fontsize='small')
            self.ax4.plot(self.wave, [
                self.data[i][self.origin[1]][self.origin[0]]
                for i in range(self.data.shape[0])
            ],
                          color=color,
                          linewidth=0.8)
            self.key2 = self.fig.canvas.mpl_connect("key_press_event",
                                                    self.keyPress2)
            plt.draw()
            self.ax4.set_ylim(self.ymin, self.ymax)
        elif event.key == 'p':
            self.origin = (round(event.xdata), round(event.ydata))
            print('Spectrum has been saved using FWHM ch3 PSF extraction area', event.key)
            print('x = %i, y = %i, PSF extract' % (self.origin[0], self.origin[1]))
            fwhm_psf3 = 0.5 # FWHM of PSF in channel 3
            n_pixels = np.round(fwhm_psf3/self.pixscale)
            color = [random.random(), random.random(), random.random()]
            square = plt.Rectangle(
                (self.origin[0] - n_pixels, self.origin[1] - n_pixels), n_pixels*2, n_pixels*2, fc=color, alpha=0.7)
            self.ax1.add_patch(square)
            self.ax1.text(self.origin[0]+1, self.origin[1]+1, s="(%i, %i)"%(self.origin[0], self.origin[1]), fontsize='small')
            flux = self.data[:,int(self.origin[1]-n_pixels):int(self.origin[1]+n_pixels),int(self.origin[0]-n_pixels):int(self.origin[0]+n_pixels)]*self.pixscale**2
            if 'sr' in self.unit: 
                flux /= 4.255e10 # flux in MJy
                flux /= 1e6 # flux in Jy
            else: 
                flux *= 1e-6 # flux in Jy
            mean_flux = np.zeros(flux.shape[0])
            for i in range(flux.shape[0]):
                mean_flux[i] = np.mean(flux[i])
            self.ax4.plot(self.wave,
                mean_flux,
                color=color,
                linewidth=0.8)
            self.key2 = self.fig.canvas.mpl_connect("key_press_event",
                                                    self.keyPress2)
            plt.draw()
            # self.ax4.set_ylim(self.ymin, self.ymax)
        

    def keyPress2(self, event):
        '''You can assign here different keys to run functions inside the lower right panel
        '''
        # sys.stdout.flush()
        if event.key == '<':
            ix = event.xdata
            print('Create a function to execute -> e.g. PAHFIT')

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        explorer = Explorer(filename=sys.argv[1])
    else:
        print('No cube found. Try command: python explore_cube.py path_cube.fits')
    explorer.start()
