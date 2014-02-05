#!/usr/bin/perl

package mapping;

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
    ServiceName => 'mapping', # any name you like  
    Description => 'mapping', 
    InputClass => 'file:///media/marco/Elements/sadi/mapping.owl#file', # defines a single service input
    OutputClass => 'file:///media/marco/Elements/sadi/mapping.owl#bam_file', # defines a single service output
    Authority => 'sadiframework.org', # domain of organization providing service
    Provider => 'myaddress@organization.org', # contact e-mail address of service provider
};

my $service = mapping->new(%$config);
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

    my $fastq_1_property = RDF::Trine::Node::Resource->new('file:///media/marco/Elements/sadi/mapping.owl#fastq_1');
    my $fastq_2_property = RDF::Trine::Node::Resource->new('file:///media/marco/Elements/sadi/mapping.owl#fastq_2');
    my $adapter_property = RDF::Trine::Node::Resource->new('file:///media/marco/Elements/sadi/mapping.owl#adapter');
    my $genome_property = RDF::Trine::Node::Resource->new('file:///media/marco/Elements/sadi/mapping.owl#genome');
    my $offset_property = RDF::Trine::Node::Resource->new('file:///media/marco/Elements/sadi/mapping.owl#offset');
    my $mapping_property = RDF::Trine::Node::Resource->new('file:///media/marco/Elements/sadi/mapping.owl#mapping');

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

    my ($fastq_1) = $input_model->objects( @$inputs[0],$fastq_1_property);  
    my ($fastq_2) = $input_model->objects( @$inputs[0],$fastq_2_property);  
    my ($adapter) = $input_model->objects( @$inputs[0],$adapter_property);
    my ($genome) = $input_model->objects( @$inputs[0],$genome_property);
    my ($offset) = $input_model->objects( @$inputs[0],$offset_property);
    #`fastq-mcf -o $fastq_1"_trimmed" -o $fastq_2"_trimmed" -0 -l 17 -u $adapter $fastq_1 $fastq_2  `;
    #`gmap_build -d genome -D ./genome $genome`;
    my $sam_file = $fastq_1->value;
    $sam_file =~ s/_1//;
    $sam_file =~ s/fastq/sam/;
    #`gsnap -B 5 -t 4 -A sam -d genome -D ./genome $fastq_1"_trimmed" $fastq_2"_trimmed"  > $sam_file`;
    my $bam_file_1 = $sam_file;
    $bam_file_1 =~ s/.sam/_1.bam/;
    my $bam_sorted = $sam_file;
    $bam_sorted =~ s/.sam/.sorted/;
    #`samtools view -bS -h -f 0x0040 $sam_file > $bam_file_1`;
    #`samtools sort $bam_file_1 $bam_sorted`;
    #`samtools index $bam_sorted".bam"`;
    my $filt_file = $sam_file;
    $filt_file =~ s/.sam//;
    #`python filter.py $bam_sorted".bam" Magnaporthe_oryzae.MG8.18.dna.toplevel.fa $offset $filt_file`;
    my $bam_file = $filt_file . ".bam";

    my $mapping = sprintf("Final alignment will be available in '%s'", $bam_file);
    my $mapping_literal = RDF::Trine::Node::Literal->new($mapping);
    my $statement = RDF::Trine::Statement->new(@$inputs[0], $mapping_property, $mapping_literal);
    $output_model->add_statement($statement);
    
}    


