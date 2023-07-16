import vcfpy
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i')
parser.add_argument('--output', '-o')
args = parser.parse_args()

reader = vcfpy.Reader.from_path(args.input)
writer = vcfpy.Writer.from_path(args.output, reader.header)

for rec in reader:
    ref = rec.REF
    alt = rec.ALT[0].value
    if len(ref) == len(alt):
        writer.write_record(rec)
    else:
        rf, rw, af, aw = rec.INFO['SR']
        if rf + rw < af + aw:
            writer.write_record(rec)