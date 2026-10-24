function LocusID = protein_name_to_bsu(ProteinName,translation_table)

LocusID = repmat({''},length(ProteinName),1);

% add hypothetical BSU ID for protein complex SigK = spoIIIC + SpoIVCB
%translation_table = [translation_table; {'sigK','BSU25760_BSU25760','sigK','sigK'}];

translation_BSU   = translation_table(:,2);
translation_name  = translation_table(:,1);
ll                = label_names(lower(ProteinName),lower(translation_name));
LocusID(find(ll)) = translation_BSU(ll(find(ll)));
