import textgrid

FileName = "..\\new\\ydy\\Language\\6_13\\6_ansi.TextGrid"
# We must save the test.TextGrid file in ANSI and rename it before open it

with open(FileName, mode = "r") as FileIn:
    Raw = textgrid.TextGrid(FileIn.read())

# Raw.tiers is a list of Tiers
# For example, Tier No.2 is Raw.tiers[1]

TierWord = Raw.tiers[1].make_simple_transcript()
TierPrpt = Raw.tiers[2].make_simple_transcript()

# The transcript of a Tier is a list of Intervals
# An Interval is a list of two numbers (the starting and ending time)
# and its text

# Here listing the Intervals in time order

pIndex = 0
strOutput = ''  # strOutput is the formatted output
for Interval in TierWord:
   if len(Interval[2]) == 0:
      continue
   startT, endT = round(float(Interval[0]), 3), round(float(Interval[1]), 3)
   strApp = '{}\t{}\t{}'.format(startT, endT, Interval[2])
   strOutput += strApp
   while pIndex < len(TierPrpt) and float(TierPrpt[pIndex][0]) < float(Interval[0]):
      pIndex = pIndex + 1
   while pIndex < len(TierPrpt) and float(TierPrpt[pIndex][1]) <= float(Interval[1]):
      strApp = '\t{}'.format(TierPrpt[pIndex][2])
      strOutput += strApp
      pIndex = pIndex + 1
   strOutput += '\n'

fout = open('..\\new\\ydy\\Language\\6_13\\Label.txt', 'w')
fout.write(strOutput)
fout.close()
