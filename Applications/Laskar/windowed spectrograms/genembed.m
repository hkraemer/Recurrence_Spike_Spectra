function YY = genembed(Y, taus, ts)
    assert(taus(1) == 0, "provide a delay-vector, which has a 0 in its first entry.")
    assert(length(taus)==length(ts), "provide delay- and ts-index vectors, which have the same length.")
    if size(Y,2)>size(Y,1)
        Y = Y';
    end
    assert(max(ts) <= size(Y,2), "The maximum ts-index must be smaller or equal to the number of provided time series.")
    
    YY = Y(:,ts(1));
    for i = 2:length(taus)
        YY = embed_for_prediction(YY, Y(:,ts(i)), taus(i));
    end
end

function Y2 = embed_for_prediction(Y, x, tau)
    N = length(Y);
    MM = length(x);
    MMM = MM - tau;
    M = min([N, MMM]);
    Y2 = horzcat(Y(end-M+1:end,:), x(end-M-tau+1:end-tau));
end