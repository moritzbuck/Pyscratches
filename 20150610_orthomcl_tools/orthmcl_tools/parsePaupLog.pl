#!/usr/bin/perl -w

use strict;

my ($pauplog)  = @ARGV;
my $nodesfile  = $pauplog .'.nodes';
my $treefile  = $pauplog .'.tree';
my $changefile = $pauplog .'.changes';


# Get the true values from the symbols
my @symbols = (0..9,'A'..'I');
my %STATES;
for(my $i=0;$i<@symbols;$i++){
    $STATES{$symbols[$i]}=$i;
#    print "$symbols[$i] $STATES{$symbols[$i]}\n";
}


open(IN, "$pauplog");
# Read away start
while(<IN>){ last if /Data matrix and reconstructed states/;}
<IN>;<IN>;<IN>;

####################################
# Get matrix including all nodes
open(NODES, ">$nodesfile");
open(NODES2, ">$nodesfile.johan");
my %SEQS;
my %SYMSEQS;
while(<IN>){
    if(/^\-\-/){
	while(<IN>){
	    last if /^\s*$/;
	    my ($seqid, $seq) =split();
	    my (@chars) = split(//,$seq);
	    my @copynumber;
	    foreach my $char (@chars){
		push(@copynumber, $STATES{$char});
	    }  
	    push( @{$SEQS{$seqid}}, @copynumber);
	    push( @{$SYMSEQS{$seqid}}, @chars);
	}
    }
    last if /Tree length/;
}

foreach my $seqid (sort keys %SEQS){
    print NODES ">$seqid\n";
    print NODES join(" ",@{$SEQS{$seqid}}), "\n";
}

foreach my $seqid (sort keys %SYMSEQS){
    print NODES2 "$seqid  ";
    print NODES2 join("",@{$SYMSEQS{$seqid}}), "\n";
}
close(NODES);
close(NODES2);
####################################


####################################
# Print the tree
open(TREE, ">$treefile");
while(<IN>){
    last if /Character change list/;
    print TREE;
}
close(TREE);
####################################
# Read away junk
while(<IN>){last if /^\-\-/;}



####################################w
# Get the changes

open(CHANGE, ">$changefile");
print CHANGE "#branch\torthoMCL\tfromstate\ttostate\tdiff\tprecise\n";
my $current_faMCL;
my %SPLITS;
while(<IN>){
    last if /^\s*$/;
    my ($penalty,$fromnode,$fromstate,$arrow,$tostate,$tonode);
    if(/^[0-9]+/){
	(my $clust,my $CI,$penalty,$fromnode,$fromstate,$arrow,$tostate,$tonode) = split();
	$current_faMCL = '' . sprintf("%d",$clust);
    }else{
	($penalty,$fromnode,$fromstate,$arrow,$tostate,$tonode) = split();
    }
    my $splitid = $fromnode . '_to_' . $tonode;
    my $diff = $STATES{$tostate} - $STATES{$fromstate};
    my $precise = 'precise';
    $precise = 'nonprecise' if ($fromstate eq $symbols[-1] || $tostate eq $symbols[-1]);
    $SPLITS{$splitid}->{$current_faMCL}->{from} = $STATES{$fromstate};
    $SPLITS{$splitid}->{$current_faMCL}->{to} = $STATES{$tostate};
    $SPLITS{$splitid}->{$current_faMCL}->{diff} = $diff;
    $SPLITS{$splitid}->{$current_faMCL}->{precise} = $precise;

}

foreach my $splitid (sort keys %SPLITS){
    foreach my $faMCL (sort keys %{$SPLITS{$splitid}}){
	print CHANGE "$splitid\t$faMCL\t";
	print CHANGE $SPLITS{$splitid}->{$faMCL}->{from}, "\t";
	print CHANGE $SPLITS{$splitid}->{$faMCL}->{to}, "\t";
	print CHANGE $SPLITS{$splitid}->{$faMCL}->{diff}, "\t";
	print CHANGE $SPLITS{$splitid}->{$faMCL}->{precise}, "\n";
    }
}
close(CHANGE);
####################################w


    
