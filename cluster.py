import sys, numpy as np
import scipy.stats, math


def median(lst):
    return np.median(np.array(lst))

def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med))

class Cluster:
   def __init__(self):
     self.chrx = None
     self.start = None
     self.end = None
     self.height = None  
     self.expr = None
     self.ratio = None
     self.rpkm = None
     self.zscore = None 

coverage = int(sys.argv[2])
current_cluster = None
current_chrx = None
clusters = []
reading = False
for line in open(sys.argv[1]):
  chrx, pos, val = line.strip().split('\t')
  pos = int(pos)
  val = float(val)

  if reading == False:
    if val == 0:
      continue
    else:
      cluster = Cluster()
      cluster.chrx = chrx
      cluster.start = pos
      cluster.end = pos
      cluster.height = val
      cluster.expr = val
      current_cluster = cluster
      current_chrx = chrx
      clusters.append(cluster)      
      reading = True
  else:
    if chrx == current_chrx:
      if val > 0 :
        current_cluster.end = pos
        cluster.expr += val
        if val >  current_cluster.height:
           current_cluster.height = val
      else:
        reading = False
        current_cluster = None
    else:
      current_chrx = chrx
      reading = False
      current_cluster = None

#rpkms = []
#for cluster in clusters:
	#rpkms.append(((cluster.expr  * 1000000.0) / (coverage * ((cluster.end - cluster.start) / 1000.0))))
	

best_lambda = -2
best_skew = 100
for test_lambda in np.arange(-2,2,0.1):
  #ratios = []
  rpkms = []
  for cluster in clusters:
    #cluster.ratio = ((float(cluster.height)   /  (cluster.end - cluster.start)) ** test_lambda - 1) / test_lambda
    cluster.rpkm = (((cluster.expr  * 1000000.0) / (coverage * ((cluster.end - cluster.start) / 1000.0))) ** test_lambda - 1) / test_lambda
    #ratios.append( cluster.ratio )
    rpkms.append(cluster.rpkm)
    #print round(cluster.ratio, 2)
  #if abs(scipy.stats.skew(ratios)) < abs(best_skew): 
  #  best_lambda =  test_lambda
  #  best_skew = scipy.stats.skew(ratios)
  if abs(scipy.stats.skew(rpkms)) < abs(best_skew): 
    best_lambda =  test_lambda
    best_skew = scipy.stats.skew(rpkms)

#ratios = []
rpkms = []
for cluster in clusters:
  #cluster.ratio = ((float(cluster.height)   /  (cluster.end - cluster.start)) ** best_lambda - 1) / best_lambda
  cluster.rpkm = (((cluster.expr * 1000000.0) / (coverage * ((cluster.end - cluster.start) / 1000.0))) ** best_lambda - 1) / best_lambda
  rpkms.append( cluster.rpkm ) 
  #print round(cluster.ratio, 2)
  #print cluster.rpkm

#ratios_median = median(ratios)
#ratios_mad = 1.4826*mad(ratios)
#print ratios_median,ratios_mad, scipy.stats.skew(ratios)

rpkms_median = median(rpkms)
rpkms_mad = 1.4826*mad(rpkms)
#print rpkms_median,rpkms_mad, scipy.stats.skew(rpkms)

zscores = []
for cluster in clusters:
  #cluster.zscore = (cluster.ratio - ratios_median) / ratios_mad
  cluster.zscore = (cluster.rpkm - rpkms_median) / rpkms_mad
  zscores.append(cluster.zscore)
  #print cluster.chrx, cluster.start, cluster.end, cluster.expr, cluster.rpkm, cluster.zscore

#threshold = np.percentile(zscores, 99.5)
for cluster in clusters:
  if cluster.zscore > 2.5:
    print cluster.chrx, cluster.start, cluster.end, cluster.expr, cluster.rpkm, cluster.zscore

 


