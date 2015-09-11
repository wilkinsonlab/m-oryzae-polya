import sys, socket, os


if not os.path.exists(os.getcwd() + "/" + sys.argv[2]):
    os.mkdir(os.getcwd() + "/" + sys.argv[2]) 

sck = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
sck.connect(("127.0.0.1", 60151))
sck.send("new\n")
sck.recv(1024)
sck.send("genome guy11\n")
sck.recv(1024)
sck.send("load /media/marco/Elements/EXP5/GUY11/WT_2.sorted.bam\n")
sck.recv(1024)
sck.send("load /media/marco/Elements/EXP5/GUY11/exp5_2.sorted.bam\n")
sck.recv(1024)
sck.send("load /media/marco/Elements/EXP5/GUY11/rbp35_2.sorted.bam\n")
sck.recv(1024)
#sck.send("load /media/marco/Elements/3Tfill/oryzae_21/WT-CM-3_plus.bedgraph\n")
#sck.recv(1024)
#sck.send("load /media/marco/Elements/3Tfill/oryzae_21/2D4-CM-3_plus.bedgraph\n")
#sck.recv(1024)
#sck.send("load /media/marco/Elements/3Tfill/oryzae_21/WT-CM-3_minus.bedgraph\n")
#sck.recv(1024)
#sck.send("load /media/marco/Elements/3Tfill/oryzae_21/2D4-CM-3_minus.bedgraph\n")
#sck.recv(1024)
sck.send("expand Gene\n")
sck.recv(1024)
sck.send("maxPanelHeight 500\n")
sck.recv(1024)
sck.send("snapshotDirectory " + os.getcwd() + "/" + sys.argv[2] + "\n")
sck.recv(1024)
import time
#time.sleep(10)
for line in open(sys.argv[1]):
    chrx, start, end, name = line.strip().split("\t")
    left = int(start) - 1000
    right = int(end) + 1000
    sck.send("goto " + chrx + ":" + str(left) + "-" + str(right) + "\n")
    sck.recv(1024)
    sck.send("snapshot  " + name + "_" + chrx + ":" + start + "-" + end + ".jpg\n")
    sck.recv(1024)
    
    
