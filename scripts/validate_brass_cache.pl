#!/usr/bin/env perl
use strict;
use warnings;
use IO::Uncompress::Gunzip qw($GunzipError);
use JSON::PP;
use Sanger::CGP::Vagrent::Data::Transcript;
use Sanger::CGP::Vagrent::Data::Exon;
use Sanger::CGP::Vagrent::Data::GenomicRegion;
use Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource;
use Sanger::CGP::Brass::Implement;

my ($bundle, $reference) = @ARGV;
die "Usage: $0 bundle reference.fa\n" unless defined $reference;
open my $fh, '<', "$bundle/adaptation.json" or die $!;
my $report = decode_json(do { local $/; <$fh> });
close $fh;
die "Adaptation not complete" unless $report->{status} eq 'adapted_requires_native_reader_validation';
my $expected = 0;
$expected += $_ for values %{$report->{resources}->{'vagrent.cache.gz'}->{kept_rows_by_contig}};
my %length;
open $fh, '<', "$reference.fai" or die $!;
while (<$fh>) { my @f = split /\t/; $length{$f[0]} = $f[1]; }
close $fh;
my $gz = IO::Uncompress::Gunzip->new("$bundle/vagrent.cache.gz", MultiStream => 1)
    or die $GunzipError;
my ($count, %first) = (0);
while (my $line = <$gz>) {
    my @f = split /\t/, $line, 7;
    # This is the pinned publisher's Perl serialization, as used by VAGrENT.
    my $VAR1;
    eval $f[6];
    die "Cannot deserialize $f[3]: $@" if $@;
    die "Wrong transcript object" unless ref($VAR1) eq 'Sanger::CGP::Vagrent::Data::Transcript';
    die "Transcript accession mismatch" unless $VAR1->getAccession eq $f[3];
    my @exons = $VAR1->getExons;
    die "No exons: $f[3]" unless @exons;
    for my $exon (@exons) {
        die "Embedded chromosome mismatch" unless $exon->getChr eq $f[0];
        die "Exon outside reference" unless $exon->getMinPos >= 1 &&
            $exon->getMaxPos <= $length{$f[0]} && $exon->getMinPos <= $exon->getMaxPos;
    }
    $first{$f[0]} //= [$f[1]+1, $f[2], $f[3]];
    $count++;
}
die "Truncated cache: $count != $expected" unless $count == $expected;
my $reader = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new(
    cache => "$bundle/vagrent.cache.gz", sorted => 1);
for my $chr (sort keys %first) {
    my ($start, $end, $accession) = @{$first{$chr}};
    my $region = Sanger::CGP::Vagrent::Data::GenomicRegion->new(
        species => 'human', genomeVersion => 'GRCh37', chr => $chr,
        minpos => $start, maxpos => $end);
    my @transcripts = $reader->getTranscripts($region);
    my @found = grep { $_->getAccession eq $accession } @transcripts;
    die "Native indexed lookup failed for $chr $accession" unless @found;
    die "Missing transcript sequence $accession" unless length($found[0]->getcDNASeq // '') > 0;
}
my @valid = Sanger::CGP::Brass::Implement::valid_seqs({
    centtel => "$bundle/CentTelo.tsv", genome => $reference});
my @expected_chrs = (map { "chr$_" } 1..22, 'X', 'Y');
die "BRASS chromosome selection mismatch" unless join(',', @valid) eq join(',', @expected_chrs);
print encode_json({status => 'passed', deserialized_transcripts => $count,
    native_query_contigs => [sort keys %first], brass_calling_contigs => \@valid}), "\n";
