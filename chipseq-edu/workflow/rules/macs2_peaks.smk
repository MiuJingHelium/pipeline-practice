rule call_peak:
	input:
		treatment=rules.bam_sort.output.bam
	output:
		xls="results/macs2/{sample}_{genome}_peaks.xls",   
                narrowPeak="results/macs2/{sample}_{genome}_peaks.narrowPeak" 
	log:
		"logs/macs2/{sample}_{genome}_callpeak.log"
	params:
		"-f BAM -g hs"
#	benchmark:
#		"logs/benchmarks/callpeak/{sample}_{genome}.txt"
	wrapper:
		"v2.2.1/bio/macs2/callpeak"

rule callpeak_options:
    input:
        treatment=rules.bam_sort.output.bam   # required: treatment sample(s)
              # optional: control sample(s)
    output:
        # all output-files must share the same basename and only differ by it's extension
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/macs2/callpeak.html.
        multiext("results/macs2/{sample}_{genome}",
                 "_peaks.xls",   ### required
                 ### optional output files
                 # these output extensions internally set the --bdg or -B option:
                 "_treat_pileup.bdg",
                 "_control_lambda.bdg",
                 # these output extensions internally set the --broad option:
                 "_peaks.broadPeak",
                 "_peaks.gappedPeak"
                 )
    log:
        "logs/macs2/{sample}_{genome}_callpeak.log"
    params:
        "-f BAM -g hs"
#   benchmark:
#	"logs/benchmarks/callpeak/{sample}_{genome}.txt"	
    wrapper:
        "v2.2.1/bio/macs2/callpeak"
