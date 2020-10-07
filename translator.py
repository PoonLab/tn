import io 
import re
import argparse

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-i', type=str, default=(),    
    )

    parser.add_argument(
        '-o', type=str, default=(),    
    )

    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    iFile = open(args.i, "r")
    oFile = open(args.o, "w")

    pattern = "GTRRates"

    for line in iFile:
        if re.search(pattern, line):
            s = line.split()
            oFile.write("{")
            oFile.write("\n\t\"empirical_frequencies\": true,")
            oFile.write("\n\t\"datatype\": \"DNA\", ")
            oFile.write("\n\t\"subs_model\": \"GTR\",")
            oFile.write("\n\t\"program\": \"FastTree version 2.1.10 Double precision (NoSSE3)\",")
            oFile.write("\n\t\"ras_model\": \"gamma\",")
            oFile.write("\n\t\"gamma\": {")
            oFile.write("\n\t\t\"alpha\": 1.0,")
            oFile.write("\n\t\t\"n_cats\": 20")
            oFile.write("\n\t},")
            oFile.write("\n\t\"subs_rates\": {")
            oFile.write("\n\t\t\"ac\": {},".format(s[1]))
            oFile.write("\n\t\t\"gt\": {},".format(s[2]))
            oFile.write("\n\t\t\"at\": {},".format(s[3]))
            oFile.write("\n\t\t\"ag\": {},".format(s[4]))
            oFile.write("\n\t\t\"cg\": {},".format(s[5]))
            oFile.write("\n\t\t\"ct\": {}".format(s[6]))
            oFile.write("\n\t}")
            oFile.write("\n}")
