#!/usr/bin/perl

package identify_polyA;

use strict;
use warnings;

#-----------------------------------------------------------------
# CGI HANDLER PART
#-----------------------------------------------------------------

use Log::Log4perl qw(:easy);
use base 'SADI::Simple::SyncService'; # or 'SADI::Simple::SyncService'

# send log messages at level WARN or above to STDERR
Log::Log4perl->easy_init($WARN);

my $config = {
    ServiceName => 'identify_polyA', # any name you like  
    Description => 'identify_polyA', 
    InputClass => 'file:///media/marco/Elements/sadi/identify_polyA.owl#Bam_file', # defines a single service input
    OutputClass => 'file:///media/marco/Elements/sadi/identify_polyA.owl#Bedgraph', # defines a single service output
    Authority => 'sadiframework.org', # domain of organization providing service
    Provider => 'myaddress@organization.org', # contact e-mail address of service provider
};

my $service = identify_polyA->new(%$config);
$service->handle_cgi_request;

#-----------------------------------------------------------------
# SERVICE IMPLEMENTATION PART
#-----------------------------------------------------------------

use RDF::Trine::Node::Resource;
use RDF::Trine::Node::Literal;
use RDF::Trine::Statement;

=head2 process_it

 Function: implements the business logic of a SADI service
 Args    : $inputs - ref to an array of RDF::Trine::Node::Resource
           $input_model - an RDF::Trine::Model containing the input RDF data
           $output_model - an RDF::Trine::Model containing the output RDF data
 Returns : nothing (service output is stored in $output_model)

=cut

sub process_it {

    my ($self, $inputs, $input_model, $output_model) = @_;

    my $bam_file_property = RDF::Trine::Node::Resource->new('file:///media/marco/Elements/sadi/identify_polyA.owl#bam_file');
    my $annotation_file_property = RDF::Trine::Node::Resource->new('file:///media/marco/Elements/sadi/identify_polyA.owl#annotation_file');
    my $offset_property = RDF::Trine::Node::Resource->new('file:///media/marco/Elements/sadi/identify_polyA.owl#offset');
    my $identify_polyA_property = RDF::Trine::Node::Resource->new('file:///media/marco/Elements/sadi/identify_polyA.owl#identify_polyA');

    my ($i) = 0;
    my (@fastq_files) = ("", "") ;
    foreach my $input (@$inputs) {
       
        #Log4perl 'easy mode' routines: TRACE, DEBUG, INFO, WARN, ERROR, FATAL, ALWAYS
        #INFO(sprintf('processing input %s', $input->uri));
        #my ($name) = $input_model->objects($input, $name_property);
        
        #if (!$name || !$name->is_literal) {
        #    WARN(sprintf('skipping input %s, doesn\'t have a <%s> property with a literal value', $input->uri, $name_property->uri));
        #    next;
        #}
        #my $mapping = sprintf("Aligning '%s' ...", $name->value);
        #my $mapping_literal = RDF::Trine::Node::Literal->new($mapping);

        #my $statement = RDF::Trine::Statement->new($input, $mapping_property, $mapping_literal);
        #$output_model->add_statement($statement);
    }

    my ($bam_file) = $input_model->objects( @$inputs[0],$bam_file_property);  
    my ($annotation_file) = $input_model->objects( @$inputs[0],$annotation_file_property);
    my ($offset) = $input_model->objects( @$inputs[0],$offset_property);
    my $assign_file = $bam_file->value;
    my $notassign_file = $bam_file->value;
    $assign_file =~ s/.bam/.assign/;
    $notassign_file =~ s/.bam/.notassign/;
    `python assign.py $annotation_file $bam_file $assign_file $notassign_file $offset`;
    
    my $polyA_file = $bam_file->value;
    my $notpolyA_file = $bam_file->value;
    $polyA_file =~ s/.bam/.polyA/;
    $notpolyA_file =~ s/.bam/.notpolyA/;
    system(" cat $assign_file | cut -f 2,3,4,5,6,7,8 | awk '{  if ( " . '$4' . " == \"+\" ) " . 'print $2,$1,$4,$5,$6,$7; else print $3,$1,$4,$5,$6,$7' . "  }' | sort -k4 | uniq -c > $polyA_file ");
    system(" cat $notassign_file | cut -f 2,3,4,5,6,7,8 | awk '{  if ( " . '$4' . " == \"+\" ) " . 'print $2,$1,$4,$5,$6,$7; else print $3,$1,$4,$5,$6,$7'  . "}' | sort -k4 | uniq -c > $notpolyA_file ");
    
    my $expr_file = $bam_file->value;
    $expr_file =~ s/.bam/.expr/;
    system(" cat Magnaporthe_oryzae.MG8.18.gff3 | awk '{if(" . '$3' . " == \"gene\") print " . '$0' . "}' | grep  \"ID=.*;\"  -o | sed -e 's/ID=//' -e 's/;//' > _t");
    system(" cut -f 6 $assign_file | cat - _t | sort | uniq -c | awk '{print " . '$2"\t"$1-1' . "}' > $expr_file  ");
    system("rm _t"); 
    
    my $polyA_all_file = $polyA_file;
    my $notpolyA_all_high_file = $polyA_file;
    my $notpolyA_all_low_file = $polyA_file;
    $polyA_all_file =~ s/.polyA/.polyA_all_m/;
    $notpolyA_all_high_file =~ s/.polyA/.notpolyA_all_m_high/;
    $notpolyA_all_low_file =~ s/.polyA/.notpolyA_all_m_low/;
    system(" python polyA_extract.py $annotation_file $polyA_file $expr_file 33 0.05 all > $polyA_all_file ");
    system(" python polyA_extract_not.py $notpolyA_file 33 100 > $notpolyA_all_high_file "); 
    system(" python polyA_extract_not.py $notpolyA_file 33 10 > $notpolyA_all_low_file "); 
    

    my $identify_polyA = sprintf("PolyA identified available in '%s', '%s', '%s', '%s'", $polyA_all_file, $notpolyA_all_high_file, $notpolyA_all_low_file );
    my $identify_polyA_literal = RDF::Trine::Node::Literal->new($identify_polyA);
    my $statement = RDF::Trine::Statement->new(@$inputs[0], $identify_polyA_property, $identify_polyA_literal);
    $output_model->add_statement($statement);
    
}    
