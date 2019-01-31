import sys
from reluplexFunctions import checkRegionConfig

if __name__== "__main__":
    cfg_file = "./reluplex.cfg"
    if len(sys.argv)>1:
        cfg_file = sys.argv[1]
    checkRegionConfig(cfg_file)