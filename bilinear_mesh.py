from pylab import *

def bilinear_interpolate(im, x, y):

    """
    Taken from
    https://stackoverflow.com/questions/12729228/simple-efficient-bilinear-interpolation-of-images-in-numpy-and-python
    """
    
    x = np.asarray(x)
    y = np.asarray(y)

    x0 = np.floor(x).astype(int)
    x1 = x0 + 1
    y0 = np.floor(y).astype(int)
    y1 = y0 + 1

    x0 = np.clip(x0, 0, im.shape[1]-1);
    x1 = np.clip(x1, 0, im.shape[1]-1);
    y0 = np.clip(y0, 0, im.shape[0]-1);
    y1 = np.clip(y1, 0, im.shape[0]-1);

    Ia = im[ y0, x0 ]
    Ib = im[ y1, x0 ]
    Ic = im[ y0, x1 ]
    Id = im[ y1, x1 ]

    Ia[y0==im.shape[0]-1] = np.nan
    Ia[y0==0] = np.nan

    wa = (x1-x) * (y1-y)
    wb = (x1-x) * (y-y0)
    wc = (x-x0) * (y1-y)
    wd = (x-x0) * (y-y0)

    return wa*Ia + wb*Ib + wc*Ic + wd*Id

def get_polar(new_rmin, new_rmax, new_tmin, new_tmax,
              rcoords, tcoords, field, pxx=1000, pxy=1000):
    
    new_r = np.linspace(new_rmin,new_rmax,pxy)
    new_t = np.linspace(new_tmin,new_tmax,pxx)

    NP,NR = meshgrid(new_t, new_r)

    # Aussming planet is at rp=1 phip=0
    TM = arctan2( NR*sin(NP),(NR*cos(NP) + 1.0))
    RM = sqrt(NR**2+1.0+ 2*NR*cos(NP) )
    
    IR = interp(RM,rcoords,arange(rcoords.shape[0]))
    IT = interp(TM,tcoords,arange(tcoords.shape[0]))

    extent = [new_t.min(), new_t.max(), new_r.max(), new_r.min()]
    return  bilinear_interpolate(field, IT.ravel(), IR.ravel()).reshape(pxy,pxx), extent, new_t, new_r
