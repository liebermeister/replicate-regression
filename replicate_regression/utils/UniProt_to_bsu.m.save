function LocusID = UniProt_to_bsu(UniProtnumber,translation_table)

LocusID = repmat({''},length(UniProtnumber),1);

translation_BSU   = translation_table(:,2);
translation_name  = translation_table(:,3);
ll                = label_names(lower(UniProtnumber),lower(translation_name));
LocusID(find(ll)) = translation_BSU(ll(find(ll)));
