# usage: awk -f changeChrNamesInMODFile.awk in.mod > out.mod
{IFS="\t"; OFS="\t";
if (/^#/)
	print $0;
else
	print $1,"chr"$2,$3,$4;
}
