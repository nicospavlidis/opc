function L = fixLabels(labels)
%Enforces cluster labels to be in the range 1:K
%L = FIXLABELS(LABELS)

l = unique(labels);	
L = zeros(length(labels),1);
for i=1:length(l),
	L( labels==l(i) ) = i;
end
