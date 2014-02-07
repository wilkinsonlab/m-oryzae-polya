#!/usr/bin/perl

package create_bedgraph;

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
    ServiceName => 'create_bedgraph', # any name you like  
    Description => 'create_bedgraph', 
    InputClass => 'file:///media/marco/Elements/sadi/create_bedgraph.owl#Bam_file', # defines a single service input
    OutputClass => 'file:///media/marco/Elements/sadi/create_bedgraph.owl#Bedgraph', # defines a single service output
    Authority => 'sadiframework.org', # domain of organization providing service
    Provider => 'myaddress@organization.org', # contact e-mail address of service provider
};

my $service = create_bedgraph->new(%$config);
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

    my $bam_file_property = RDF::Trine::Node::Resource->new('file:///media/marco/Elements/sadi/create_bedgraph.owl#bam_file');
    my $annotation_file_property = RDF::Trine::Node::Resource->new('file:///media/marco/Elements/sadi/create_bedgraph.owl#annotation_file');
    my $offset_property = RDF::Trine::Node::Resource->new('file:///media/marco/Elements/sadi/create_bedgraph.owl#offset');
    my $create_bedgraph_property = RDF::Trine::Node::Resource->new('file:///media/marco/Elements/sadi/create_bedgraph.owl#create_bedgraph');

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
    print $assign_file, $notassign_file, $annotation_file, $bam_file;
    `python assign.py $annotation_file $bam_file $assign_file $notassign_file $offset`;
    my $bedgraph_plus_file = $bam_file->value;
    my $bedgraph_minus_file = $bam_file->value;
    $bedgraph_plus_file =~ s/.bam/_plus.bedgraph/;
    $bedgraph_minus_file =~ s/.bam/_minus.bedgraph/;
    my $bedgraph_not_plus_file = $bam_file->value;
    my $bedgraph_not_minus_file = $bam_file->value;
    $bedgraph_not_plus_file =~ s/.bam/_not_plus.bedgraph/;
    $bedgraph_not_minus_file =~ s/.bam/_not_minus.bedgraph/;
    
    system(" cat $assign_file | cut -f 2,3,4,5 | awk '{  if ( " . '$4' . " == \"-\" ) print " . '$3,$1' . " }' | sort -n | uniq -c | awk '{ printf " . '"%s\t%d\t%d\t%d\n",$3,$2-1,$2,$1' .  "}' > $bedgraph_plus_file ");
    system(" cat $assign_file | cut -f 2,3,4,5 | awk '{  if ( " . '$4' . " == \"+\" ) print " . '$2,$1' . " }' | sort -n | uniq -c | awk '{ printf " . '"%s\t%d\t%d\t%d\n",$3,$2-1,$2,$1' .  "}' > $bedgraph_minus_file ");
    system(" cat $notassign_file | cut -f 2,3,4,5 | awk '{  if ( " . '$4' . " == \"-\" ) print " . '$3,$1' . " }' | sort -n | uniq -c | awk '{ printf " . '"%s\t%d\t%d\t%d\n",$3,$2-1,$2,$1' .  "}' > $bedgraph_not_plus_file ");
    system(" cat $notassign_file | cut -f 2,3,4,5 | awk '{  if ( " . '$4' . " == \"+\" ) print " . '$2,$1' . " }' | sort -n | uniq -c | awk '{ printf " . '"%s\t%d\t%d\t%d\n",$3,$2-1,$2,$1' .  "}' > $bedgraph_not_minus_file ");


    my $create_bedgraph = sprintf("Final bedgraphs will be available in '%s', '%s', '%s', '%s'", $bedgraph_plus_file, $bedgraph_minus_file, $bedgraph_not_plus_file, $bedgraph_not_minus_file);
    my $create_bedgraph_literal = RDF::Trine::Node::Literal->new($create_bedgraph);
    my $statement = RDF::Trine::Statement->new(@$inputs[0], $create_bedgraph_property, $create_bedgraph_literal);
    $output_model->add_statement($statement);
    
}    
