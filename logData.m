function logData(src,evt,fid)

    data = [evt.TimeStamps,evt.Data]';
    fwrite(fid,data,'double');
end