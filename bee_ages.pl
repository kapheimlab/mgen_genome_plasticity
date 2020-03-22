# Author: Robert Waterhouse
#!/usr/local/bin/perl

use strict;
use REST::Client;
use JSON;
use Data::Dumper;


my $debug=1; # set here $debug=0/1 to turn on/off debug data printing
if($debug==0) { open(DEBUG,">debug_data.txt") || die $!; }

# Set connection to website:
#my $host = 'http://ezmeta/odb/dev';
my $host = 'http://www.orthodb.org/v9/';
my $client = REST::Client->new(host => $host);

# Choose levels:
my %levels = (
    'Apoidea' => 34735,
    'Aculeata' => 7434,
    'Hymenoptera' => 7399,
    'Holometabola' => 33392,
    'Insecta' => 50557,
    'Arthropoda' => 6656,
    'Metazoa' => 33208,
    'Vertebrata' => 7742
    );
undef my %level2species;
undef my %level2speciescount;
undef my %taxid2name;

# Loop through each level defined above to define species lists
foreach my $level (sort keys %levels) {
    
    $client->GET("search?limit=1&level=$levels{$level}&query=&universal=1");   
    if($client->responseCode() !=200) { print "Bad query\n"; exit(); }

    # Convert JSON to PERL data structure
    my $result = decode_json( $client->responseContent() );
    
    # for DEBUG: print JSON and/or PERL result:
    if($debug==0) { 
	print DEBUG "RESULT1: JSON\n\n";
	print DEBUG $client->responseContent() . "\n\n";
	print DEBUG "RESULT1: PERL\n\n";
	print DEBUG Dumper($result) . "\n\n";
    }

    my @groups = @{ $result->{'data'} };

    foreach my $group (@groups) {
	$client->GET("orthologs?id=$group");
	if($client->responseCode() !=200) { print "Bad query\n"; }
	my $result = decode_json( $client->responseContent() );
	
	if($debug==0) {
	    print DEBUG "RESULT2: JSON\n\n";
	    print DEBUG $client->responseContent() . "\n\n";
	    print DEBUG "RESULT2: PERL\n\n";
	    print DEBUG Dumper($result) . "\n\n";
	}

	# Loop through each species in the group
	$level2speciescount{$level}=0;
	foreach my $species ( @{ $result->{'data'} } ) {
	    $level2speciescount{$level}++;
	    push(@{$level2species{$level}},$species->{'organism'}{'id'});
	    $taxid2name{$species->{'organism'}{'id'}}=$species->{'organism'}{'name'};
	}
    }
    
}

$taxid2name{'115081'}='Megalopta genalis';
$taxid2name{'178049'}='Nomia melanderi';
foreach my $level (keys %level2species) {
    if($level eq 'Vertebrata') { next; }
    push(@{$level2species{$level}},'115081','178049');
    $level2speciescount{$level}++;
    $level2speciescount{$level}++;
}

my @levelorder=sort( { $level2speciescount{$b} <=> $level2speciescount{$a} } keys %level2speciescount);
my @tmp;
foreach my $level (@levelorder) { if($level ne 'Vertebrata') { push(@tmp,$level); } }
@levelorder=('Vertebrata',@tmp);

undef my %level2speciesUNI;
for(my $level=0; $level<$#levelorder; $level++) {
    my $levelname=$levelorder[$level];
    my $nextlevelname=$levelorder[$level+1];
    undef my %excludes;
    foreach my $species (@{$level2species{$nextlevelname}}) { $excludes{$species}=1; }
    if($levelname eq 'Metazoa') {
	foreach my $species (@{$level2species{'Vertebrata'}}) { $excludes{$species}=1; }
    }
    foreach my $species (@{$level2species{$levelname}}) {
	if(!defined($excludes{$species})) { push(@{$level2speciesUNI{$levelname}},$species); }
    }
    if($levelname eq 'Vertebrata') { @{$level2speciesUNI{$levelname}}=@{$level2species{$levelname}}; }
    if($debug==0) {
	print DEBUG "$levelname\t$levels{$levelname}\t" . join(";",sort(@{$level2species{$levelname}})) . "\n";
	print DEBUG "$levelname\t$levels{$levelname}\t" . join(";",sort(@{$level2speciesUNI{$levelname}})) . "\n";
    }
    if($nextlevelname eq $levelorder[$#levelorder]) {
	foreach my $species (@{$level2species{$nextlevelname}}) {
	    push(@{$level2speciesUNI{$nextlevelname}},$species);
	}
	if($debug==0) {
	    print DEBUG "$nextlevelname\t$levels{$nextlevelname}\t" . join(";",sort(@{$level2species{$nextlevelname}})) . "\n";
	    print DEBUG "$nextlevelname\t$levels{$nextlevelname}\t" . join(";",sort(@{$level2speciesUNI{$nextlevelname}})) . "\n";
	    print DEBUG "\n\n\n";
	}
    }
}


# read in the two mapped species first
# and in so doing, actually get all root-level OGs

undef my %og2genes;
undef my %gene2ogN;
open(IN,"node_33208_taxid_115081.og") || die $!;
my @lines=<IN>;
close(IN);
foreach my $line (@lines) {
    if($line=~/^#/) { next; }
    my @bits=split(/\s+/,$line);
    push(@{$og2genes{$bits[0]}},$bits[1]);
    $gene2ogN{$bits[1]}=$bits[0];
}
open(IN,"node_33208_taxid_178049.og") || die $!;
@lines=<IN>;
close(IN);
foreach my $line (@lines) {
    if($line=~/^#/) { next; }
    my @bits=split(/\s+/,$line);
    if($bits[1]=~/^178049:/) { push(@{$og2genes{$bits[0]}},$bits[1]); }
    if(defined($gene2ogN{$bits[1]})) {
	if($gene2ogN{$bits[1]}!=$bits[0]) { print "OG mismatch: $gene2ogN{$bits[1]} != $bits[0]\n"; }
    }
    else { $gene2ogN{$bits[1]}=$bits[0]; }
}
undef my %int2pub;
open(IN,"115081.fs.maptxt") || die $!;
@lines=<IN>;
close(IN);
foreach my $line (@lines) {
 chomp($line);
 my @bits=split(/\s+/,$line);
 $int2pub{$bits[0]}=$bits[1];
}
open(IN,"178049.fs.maptxt") || die $!;
@lines=<IN>;
close(IN);
foreach my $line (@lines) {
 chomp($line);
 my @bits=split(/\s+/,$line);
 $int2pub{$bits[0]}=$bits[1];
}

# get pubids for Bimpa
open(IN,"Bimpa_ODB9_pub_gene_IDmap.txt") || die $!;
@lines=<IN>;
close(IN);
foreach my $line (@lines) {
    chomp($line);
    my ($int,$pub)=split(/\s+/,$line);
    my ($a,$b,$c)=split(/\|/,$pub);
    $int2pub{$int}=$c;
}
# get proper pubids for some bees
open(IN,"bee_ids.txt") || die $!;
@lines=<IN>;
close(IN);
undef my %beefix;
foreach my $line (@lines) {
    chomp($line);
    my ($a,$b)=split(/\s+/,$line);
    $beefix{$a}=$b;
}


# now get all Metazoa groups and assess phyletic profiles

#$client->GET("search?limit=500&level=33208&query=");   
$client->GET("search?limit=100000&level=33208&query=");   
if($client->responseCode() !=200) { print "Bad query\n"; exit(); }

# Convert JSON to PERL data structure
my $result = decode_json( $client->responseContent() );

# for DEBUG: print JSON and/or PERL result:
if($debug==0) { 
    print DEBUG "RESULT3: JSON\n\n";
    print DEBUG $client->responseContent() . "\n\n";
    print DEBUG "RESULT3: PERL\n\n";
    print DEBUG Dumper($result) . "\n\n";
}

my @groups = @{ $result->{'data'} };

undef my %ogN2eog;
undef my %eog2species;
undef my %eog2description;
undef my %eog2evorate;

foreach my $group (@groups) {
    $client->GET("orthologs?id=$group");
    if($client->responseCode() !=200) { print "Bad query\n"; }
    my $result = decode_json( $client->responseContent() );
    
    if($debug==0) {
	print DEBUG "RESULT4: JSON\n\n";
	print DEBUG $client->responseContent() . "\n\n";
	print DEBUG "RESULT4: PERL\n\n";
	print DEBUG Dumper($result) . "\n\n";
    }
    
    # Loop through each species in the group
    foreach my $species ( @{ $result->{'data'} } ) {
	my $taxid=$species->{'organism'}{'id'};
	push(@{$eog2species{$group}},$taxid);
	# Loop through each gene from the species
	foreach my $gene ( @{ $species->{'genes'} } ) {
	    my $odbid=$gene->{'gene_id'}{'param'};
	    if(!defined($int2pub{$odbid})) { $int2pub{$odbid}=$gene->{'gene_id'}{'id'}; }
	    if(!defined($ogN2eog{$gene2ogN{$odbid}})) { $ogN2eog{$gene2ogN{$odbid}}=$group; }
	    else { if($ogN2eog{$gene2ogN{$odbid}} ne $group) { print "OG misclassification\n"; } }
	}
    }

    $client->GET("group?id=$group");
    if($client->responseCode() !=200) { print "Bad query\n"; }
    $result = decode_json( $client->responseContent() );
    
    if($debug==0) {
	print DEBUG "RESULT5: JSON\n\n";
	print DEBUG $client->responseContent() . "\n\n";
	print DEBUG "RESULT5: PERL\n\n";
	print DEBUG Dumper($result) . "\n\n";
    }
    
    $eog2description{$group}=$result->{'data'}{'name'};
    $eog2evorate{$group}=$result->{'data'}{'evolutionary_rate'};
    
}


open(OUT1,">groups_aged.txt") || die $!;
open(OUT2,">genes_aged.txt") || die $!;
print OUT1 "OrthoGroup";
foreach my $level (@levelorder) { print OUT1 "\t$level\[". scalar(@{$level2speciesUNI{$level}}) ."]"; }
print OUT1 "\tAge\tEvoRate\tDescription\n";
print OUT2 "OrthoGroup\tAge\tSpecies\tODBID\tGeneID\n";

foreach my $og (sort keys %og2genes) {

    if(!defined($ogN2eog{$og})) { next; }
    my $eog=$ogN2eog{$og};

    # add mapped species
    my $Mgen=0;
    my $Nmel=0;
    foreach my $gene (@{$og2genes{$og}}) {
	if($gene=~/^115081:/) { $Mgen++; }
	if($gene=~/^178049:/) { $Nmel++; }
    }
    if($Mgen>0) { push(@{$eog2species{$eog}},'115081'); }
    if($Nmel>0) { push(@{$eog2species{$eog}},'178049'); }
    

    my @clades=();
    my $bees=0;
    foreach my $level (@levelorder) {
	undef my %counts;
	foreach my $species (@{$eog2species{$eog}},@{$level2speciesUNI{$level}}) { $counts{$species}++; }
	my $shared=0;
	foreach my $species (@{$level2speciesUNI{$level}}) { if($counts{$species}==2) { $shared++;} }
	push(@clades,$shared);
	if($level eq 'Apoidea') { $bees=$shared; }
    }

    if($bees==0) { next; }

    # set type
    my $type='Apoidea';
    if($clades[0]>0) { $type='Vertebrata'; }
    elsif($clades[1]>0) { $type='Metazoa'; }
    elsif($clades[2]>0) { $type='Arthropoda'; }
    elsif($clades[3]>0) { $type='Insecta'; }
    elsif($clades[4]>0) { $type='Holometabola'; }
    elsif($clades[5]>0) { $type='Hymenoptera'; }
    elsif($clades[6]>0) { $type='Aculeata'; }
	    
    print OUT1 "$ogN2eog{$og}\t";
    print OUT1 join("\t",@clades);
    if(!defined($eog2description{$eog})) { $eog2description{$eog}=''; }
    print OUT1 "\t$type\t$eog2evorate{$eog}\t$eog2description{$eog}\n";

    undef my %selspec;
    foreach my $spec (@{$level2species{'Apoidea'}}) { $selspec{$spec}=1; }
    foreach my $gene (sort @{$og2genes{$og}}) {
	my ($tax,$gen)=split(/:/,$gene);
	if(!defined($selspec{$tax})) { next; }
	if(defined($beefix{$int2pub{$gene}})) {$int2pub{$gene}=$beefix{$int2pub{$gene}}; }
	print OUT2 "$eog\t$type\t$taxid2name{$tax}\t$gene\t$int2pub{$gene}\n";
    }

}



close(OUT1);
close(OUT2);
if($debug==0) { close(DEBUG); }
