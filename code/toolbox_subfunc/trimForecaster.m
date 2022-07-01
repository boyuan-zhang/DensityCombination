function ecbspf00 = trimForecaster(ecbspf00)

% ecbspf00 = ecbspf_infl_1y;

ns = size(ecbspf00,1)-1;

% find a set of forecaster ID
id = [];
for sind = 1:ns
   id = [id; ecbspf00(sind).id_hist];
end

% unique ID
id_uni = unique(id);

% table of presence of forecaster in each survey
tab_presence = zeros(ns, numel(id_uni));

for sind = 1:ns
   temp_id = ecbspf00(sind).id_hist;
   tab_presence(sind,:) = ismember(id_uni, temp_id);
end

% find inactive forecasters
id_active = ones(numel(id_uni),1);

for i = 1:numel(id_uni)
   
   counter = 0;
   pre = 1;
   
   for sind = 9:ns
      temp = tab_presence(sind, i);
      
      % counter number of consecutive zeros
      if temp == 0
         counter = counter + 1;
         
         if pre == 1
            counter = 1;
         end
         
      end
      
      % if miss 4, inactive
      if counter == 4
         id_active(i) = 0;
      end
      
      pre = temp;
   end
   
end

id_active = id_uni(logical(id_active));

% keep active forecasters
for sind = 1:ns
   temp_id = ecbspf00(sind).id_hist;
   tab_presence(sind,:) = ismember(id_uni, temp_id);
end

% update id and hist
for sind = 1:ns
   
   temp_id = ecbspf00(sind).id_hist;
   temp_hist = ecbspf00(sind).hist;
   
   ecbspf00(sind).id_hist = id_active;
   
   new_hist = nan*ones(numel(id_active), size(temp_hist,2));
   
   for i = 1:numel(id_active)
      if ismember(id_active(i), temp_id)
         
         temp_pos = temp_id == id_active(i);
         new_hist(i,:) = temp_hist(temp_pos,:);
         
      end
   end
   
   ecbspf00(sind).hist = new_hist;
   
   ecbspf00(sind).hist_avg_a = mean(temp_hist,1,'omitnan');
   
   ecbspf00(sind).hist_avg_s = mean(new_hist,1,'omitnan');
   
end

end