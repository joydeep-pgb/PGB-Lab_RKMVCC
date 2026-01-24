#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);
use Parallel::ForkManager;
use POSIX qw(strftime);

# ================= CONFIGURATION =================
my $INPUT_FOLDER  = "/run/media/joydeep/One_HDD/Sorghum_MetaDEG/Salt/FASTQ/";
my $OUTPUT_FOLDER = "/run/media/joydeep/One_HDD/Sorghum_MetaDEG/Salt/FASTQ/FastQC_Results/";
my $NUM_WORKERS   = 6;
# =================================================

# ----------------- Setup output dir -----------------
unless (-d $OUTPUT_FOLDER) {
    make_path($OUTPUT_FOLDER) or die "Cannot create output directory\n";
}

# ----------------- Logging -----------------
my $log_file = "$OUTPUT_FOLDER/fastqc_analysis.log";
open(my $LOG, ">", $log_file) or die "Cannot open log file\n";

sub log_msg {
    my ($level, $msg) = @_;
    my $time = strftime("%Y-%m-%d %H:%M:%S", localtime);
    print $LOG "$time - $level - $msg\n";
}

# ----------------- Find FASTQ files -----------------
opendir(my $DIR, $INPUT_FOLDER) or die "Cannot open input directory\n";

my @fastq_files = grep {
    /\.(fastq|fq)(\.gz)?$/
} readdir($DIR);

closedir($DIR);

if (!@fastq_files) {
    die "Error: No FASTQ files found in $INPUT_FOLDER\n";
}

print "Starting FastQC on ", scalar(@fastq_files), " files...\n";
print "Using $NUM_WORKERS parallel workers.\n";

# ----------------- Parallel execution -----------------
my $pm = Parallel::ForkManager->new($NUM_WORKERS);

my $success = 0;
my $failure = 0;
my $count   = 0;
my $total   = scalar(@fastq_files);

$pm->run_on_finish(sub {
    my ($pid, $exit_code, $ident, $signal, $core_dump, $data_ref) = @_;
    $count++;

    if ($data_ref->{status} eq "OK") {
        $success++;
        log_msg("INFO", "Successfully processed: $data_ref->{file}");
    } else {
        $failure++;
        log_msg("ERROR", "Error processing $data_ref->{file}");
    }

    printf "\rProgress: %d / %d", $count, $total;
});

foreach my $file (@fastq_files) {
    $pm->start and next;

    my $full_path = "$INPUT_FOLDER/$file";
    my $cmd = "fastqc -o $OUTPUT_FOLDER $full_path";

    my $exit = system($cmd);

    if ($exit == 0) {
        $pm->finish(0, { status => "OK", file => $file });
    } else {
        $pm->finish(0, { status => "FAIL", file => $file });
    }
}

$pm->wait_all_children;
close $LOG;

# ----------------- Final summary -----------------
print "\n\n--- Analysis Completed ---\n";
print "Successfully processed: $success\n";
print "Failures: $failure\n";
print "Logs saved to: $log_file\n";
print "Reports saved to: $OUTPUT_FOLDER\n";