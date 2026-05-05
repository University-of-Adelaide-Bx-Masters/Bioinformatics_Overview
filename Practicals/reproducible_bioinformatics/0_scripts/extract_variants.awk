#!/usr/bin/awk -f

BEGIN {
    FS = "[,\t]"
    OFS = "\t"
    print "Sample","Variant","ChrPos","Ref","Alt","GT","FILTER","FT","AF","DP"
}

# First file (CSV): build lookup
FNR == NR {
    coords[$2] = $1
    next
}

# Skip headers in VCF
/^#/ { next }

{
    key = $1 ":" $2

    if (key in coords) {

        # get sample name from filename
        n = split(FILENAME, f, "/")
        sample = f[n]

        ref = $4
        alt = $5
        filter = $7

        # --- parse INFO ---
        af = "."
        dp = "."

        ninfo = split($8, info, ";")
        for (i = 1; i <= ninfo; i++) {
            split(info[i], kv, "=")
            if (kv[1] == "AF") af = kv[2]
            if (kv[1] == "DP") dp = kv[2]
        }

        # --- parse FORMAT/sample ---
        ft = "."
        gt = "."

        split($9, fmt, ":")
        split($10, val, ":")

        for (i = 1; i <= length(fmt); i++) {
            if (fmt[i] == "GT") gt = val[i]
            if (fmt[i] == "FT") ft = val[i]
            if (fmt[i] == "DP" && dp == ".") dp = val[i]
        }

        print sample, coords[key], key, ref, alt, gt, filter, ft, af, dp
    }
}
