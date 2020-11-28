import io 
import re
import argparse
import json

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-i', type=str, default=(),    
    )

    parser.add_argument(
        '-o', type=str, default=(),    
    )
    
    parser.add_argument(
        '-datatype', type=str, default="DNA",    
    )
    
    parser.add_argument(
        '-subs_model', type=str, default="GTR",    
    )

    parser.add_argument(
        '-program', type=str, default="FastTree version 2.1.10 Double precision (NoSSE3)",    
    )

    parser.add_argument(
        '-ras_model', type=str, default="gamma",    
    )


    parser.add_argument(
        '-alpha', type=str, default="1.0",    
    )

    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    iFile = open(args.i, "r")
    oFile = open(args.o, "w")

    if args.program == "FastTree version 2.1.10 Double precision (NoSSE3)":   
        pattern = "GTRRates"
        oText = {}
        for line in iFile:
            if re.search(pattern, line):
                s = line.split()
        
                oText["empirical_frequencies"] = True
                oText["datatype"] = args.datatype
                oText["subs_model"] = args.subs_model
                oText["program"] = args.program
                oText["ras_model"] = args.ras_model
                gamma = {}
                gamma["alpha"] = 1.0
                gamma["n_cats"] = 20
                oText["gamma"] = gamma
                sub_rates = {}
                sub_rates["ac"] = float(s[1])
                sub_rates["gt"] = float(s[2])
                sub_rates["at"] = float(s[3])
                sub_rates["ag"] = float(s[4])
                sub_rates["cg"] = float(s[5])
                sub_rates["ct"] = float(s[6])
                oText["subs_rates"] = sub_rates

        json.dump(oText, oFile, indent=4)


    if args.program == "IQ-TREE multicore version 1.6.12":
        pattern="Rate parameters:"
        oText = {}
        for line in iFile:
            if re.search(pattern, line):
                s = line.split()
 
                oText["empirical_frequencies"] = True
                oText["datatype"] = args.datatype
                oText["subs_model"] = args.subs_model
                oText["program"] = args.program
                oText["ras_model"] = args.ras_model
                gamma = {}
                gamma["alpha"] = 1.0
                gamma["n_cats"] = 20
                oText["gamma"] = gamma
                sub_rates = {}
                sub_rates["ac"] = float(s[3])
                sub_rates["gt"] = float(s[13])
                sub_rates["at"] = float(s[7])
                sub_rates["ag"] = float(s[5])
                sub_rates["cg"] = float(s[9])
                sub_rates["ct"] = float(s[11])
                oText["subs_rates"] = sub_rates

        json.dump(oText, oFile, indent=4)
                