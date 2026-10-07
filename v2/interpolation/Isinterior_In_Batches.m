function in = Isinterior_In_Batches(field_polygon,x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% isinterior over many points, chunked so a complex polygon cannot exhaust
%% memory on a fine grid.
%%
%% polyshape/isinterior builds one working array sized by points times
%% vertices. That is fine for the few dozen samples this repository usually
%% tests against, but a hundred-vertex field swept down to a fraction of a
%% metre puts hundreds of thousands of points against that many vertices at
%% once, which is where it runs out. Chunking by a fixed points-times-
%% vertices budget keeps every call the same size regardless of field
%% complexity or grid resolution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = reshape(x,[],1);
y = reshape(y,[],1);
in = false(size(x));
if isempty(x); return; end
if isempty(field_polygon.Vertices); return; end

vertex_count = max(1,size(field_polygon.Vertices,1));
point_budget = max(1,floor(2e6/vertex_count));
for first = 1:point_budget:numel(x)
    last = min(numel(x),first+point_budget-1);
    in(first:last) = isinterior(field_polygon,x(first:last),y(first:last));
end

end
