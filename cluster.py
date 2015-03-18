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
     self.ratio = None
     self.zscore = None 


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
      current_cluster = cluster
      current_chrx = chrx
      clusters.append(cluster)      
      reading = True
  else:
    if chrx == current_chrx:
      if val > 0 :
        current_cluster.end = pos
        if val >  current_cluster.height:
           current_cluster.height = val
      else:
        reading = False
        current_cluster = None
    else:
      current_chrx = chrx
      reading = False
      current_cluster = None

for cluster in clusters:
    if cluster.end - cluster.start > 500:
        print cluster.chrx, cluster.start, cluster.end, cluster.height
exit()        

ratios = []
best_lambda = -3
best_skew = 100
for test_lambda in np.arange(-2,2,0.1):
 ratios = []
 for cluster in clusters:
  cluster.ratio = ((float(cluster.height)   /  (cluster.end - cluster.start)) ** test_lambda - 1) / test_lambda
  #cluster.ratio = (float(cluster.height)   /  (cluster.end - cluster.start)) 
  ratios.append( cluster.ratio )
  #print round(cluster.ratio, 3)
 if abs(scipy.stats.skew(ratios)) < abs(best_skew): 
   best_lambda =  test_lambda
   best_skew = scipy.stats.skew(ratios)

ratios = []
for cluster in clusters:
  cluster.ratio = ((float(cluster.height)   /  (cluster.end - cluster.start)) ** best_lambda - 1) / best_lambda
  ratios.append( cluster.ratio ) 

ratios_median = median(ratios)
ratios_mad = 1.4826*mad(ratios)
#print ratios_median,ratios_mad, scipy.stats.skew(ratios)

zscores = []
for cluster in clusters:
  cluster.zscore = (cluster.ratio - ratios_median) / ratios_mad
  zscores.append(cluster.zscore)

threshold = np.percentile(zscores, 99)
for cluster in clusters:
  if cluster.zscore > threshold:
    print cluster.chrx, cluster.start, cluster.end, cluster.height, cluster.ratio, cluster.zscore

 


