#!/bin/bash

awk '
	# TITLE
	/^TI  - / {
		sub(/^TI  - /, "")
		title = $0
		# Read continuation lines (indented with spaces)
		while ((getline line) > 0) {
			if (line ~ /^      /) {
				sub(/^[[:space:]]+/, "", line)
				title = title " " line
			}
			else {
				break
			}
		}
	}
	
	#PMID
	/^PMID- /{
		# Remove the "PMID- " prefix and print the first line
		sub(/^PMID- /, "")
		PMID = $0
	}
	# PMID is a single line so no need for while function
                            
	#ABSTRACT 
	/^AB  - /{
		# Remove the "AB  - " prefix and print the first line
		sub(/^AB  - /, "")
		abstract = $0
		# Read continuation lines (indented with spaces) 
		# call it nextline to avoid conflicts with the previous due to repetition 
		while ((getline nextline) > 0) {
			if (nextline ~ /^      /) {
				sub(/^[[:space:]]+/, "", nextline)
				abstract = abstract " " nextline
			}
			else {
				break
			}
		}
	}
	END {  #indicating the end of reiterations 
		print PMID "\t" title "\t" abstract     # awk recognizes \t as "\t" 
	}
' $1