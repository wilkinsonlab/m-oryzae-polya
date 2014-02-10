#!/usr/bin/perl

package create_gff;

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
    ServiceName => 'create_gff', # any name you like  
    Description => 'create_gff', 
    InputClass => 'file:///media/marco/Elements/sadi/create_gff.owl#Sample', # defines a single service input
    OutputClass => 'file:///media/marco/Elements/sadi/create_gff.owl#GFF_file', # defines a single service output
    Authority => 'sadiframework.org', # domain of organization providing service
    Provider => 'myaddress@organization.org', # contact e-mail address of service provider
};

my $service = create_gff->new(%$config);
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

    my $sample_property = RDF::Trine::Node::Resource->new('file:///media/marco/Elements/sadi/create_gff.owl#sample');
    my $create_gff_property = RDF::Trine::Node::Resource->new('file:///media/marco/Elements/sadi/create_gff.owl#create_gff');

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

    my ($sample) = $input_model->objects( @$inputs[0],$sample_property);  
    
    system("  cat $sample\"-\"1.polyA_all_m $sample\"-\"2.polyA_all_m $sample\"-\"3.polyA_all_m |  sort -k 2,7 | uniq -f 1 -c | awk '{" . 'if ($1 >= 2) print 0,$3,$4,$5,$6,$7,$8}' . "' > $sample\"-\"X.polyA_all_m ");
    system("  cat $sample\"-\"1.notpolyA_all_m_high $sample\"-\"2.notpolyA_all_m_high $sample\"-\"3.notpolyA_all_m_high |  sort -k 2,7 | uniq -f 1 -c | awk '{" . 'if ($1 >= 2) print 0,$3,$4,$5,$6,$7,$8}' . "' > $sample\"-\"X.notpolyA_all_m_high ");
    system("  cat $sample\"-\"1.notpolyA_all_m_low $sample\"-\"2.notpolyA_all_m_low $sample\"-\"3.notpolyA_all_m_low |  sort -k 2,7 | uniq -f 1 -c | awk '{" . 'if ($1 >= 2) print 0,$3,$4,$5,$6,$7,$8}' . "' > $sample\"-\"X.notpolyA_all_m_low ");

    system(" cut $sample\"-\"X.polyA_all_m  -d \" \" -f 1,2,3,4,5 | awk '{ " . 'if ($4 == "+") sense = "-"; else sense = "+"; printf "%s\tmarco\tpolyA_site\t%d\t%d\t.\t%s\t.\ttranscript=%s;value=%d\n", $3, $2, $2, sense, $5, $1  ' . "}' > $sample\"_\"polyA.gff ");
    

    my $create_gff = sprintf("GFF file available in '%s'", $sample . '_polyA.gff' );
    my $create_gff_literal = RDF::Trine::Node::Literal->new($create_gff);
    my $statement = RDF::Trine::Statement->new(@$inputs[0], $create_gff_property, $create_gff_literal);
    $output_model->add_statement($statement);
    
}    
