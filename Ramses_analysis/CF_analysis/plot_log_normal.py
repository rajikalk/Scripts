import numpy as np
import matplotlib.pyplot as plt
import scipy.stats


N_bins = 50


samples   = np.random.lognormal( mean=0.3, sigma=.4, size=100 )
bins_log10 = np.logspace( np.log10( samples.min()  ), np.log10( samples.max() ), N_bins )
counts, bin_edges, ignored = plt.hist( samples, bins_log10, histtype='stepfilled', label='histogram' )

shape, loc, scale = scipy.stats.lognorm.fit( samples, floc=0 )
x_fit       = np.linspace( samples.min(), samples.max(), 100 )
samples_fit = scipy.stats.lognorm.pdf( x_fit, shape, loc=loc, scale=scale )

# calculate length of each bin (required for scaling PDF to histogram)
bins_log_len = np.zeros( bins_log10.size )
for ii in range( counts.size):
    bins_log_len[ii] = bin_edges[ii+1]-bin_edges[ii]

# get pdf-values for same intervals as histogram
samples_fit_log = scipy.stats.lognorm.pdf( bins_log10, shape, loc=loc, scale=scale )
plt.clf()
# oplot fitted and scaled PDF into histogram
plt.plot( bins_log10, np.multiply(samples_fit_log,bins_log_len)*sum(counts), label='PDF using histogram bins', linewidth=2 )
plt.xscale( 'log' )
plt.xlim( bin_edges.min(), bin_edges.max() )
plt.legend(loc=3)
plt.show()

