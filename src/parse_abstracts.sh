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
	
	# ABSTRACT 
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
	END {  # indicating the end of reiterations 
		print title "\t" abstract 
	}
' $1
