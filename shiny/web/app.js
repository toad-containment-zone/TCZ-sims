// Cane Toad Forecast -- static Leaflet front end.
// Data comes from data/*.json, written by shiny/export_web_data.R.
(function () {
  'use strict';

  // ---- Colour bins: same as leaflet R's colorBin("YlOrRd", bins = 0,.2,...,1) ----
  // Intervals are left-closed ([0,.2), [.2,.4) ... [.8,1]), as colorBin produces.
  var BREAKS = [0, 0.2, 0.4, 0.6, 0.8, 1];
  var COLOURS = ['#FFFFB2', '#FECC5C', '#FD8D3C', '#F03B20', '#BD0026'];
  function colourFor(p) {
    for (var i = COLOURS.length - 1; i > 0; i--) if (p >= BREAKS[i]) return COLOURS[i];
    return COLOURS[0];
  }

  var GREY = { radius: 3, fillColor: '#aaaaaa', fillOpacity: 0.4 };
  var PLAY_MS = 900;

  // ---- Tabs ----
  var tabs = document.querySelectorAll('[data-tab]');
  var map; // set below; resized when the Forecast tab is shown again
  // If the map is first fitted while the Forecast tab is hidden (page opened on
  // #about), its container has zero size and the zoom comes out wrong; remember
  // the bounds and fit again once the tab is visible.
  var pendingBounds = null;
  function fitMap(bounds) {
    if (map.getContainer().clientWidth === 0) { pendingBounds = bounds; return; }
    pendingBounds = null;
    // Leaflet caches the container size and invalidateSize() is a no-op until
    // the map has a view, so give it a provisional view, refresh the size, fit.
    if (!map._loaded) map.setView(L.latLngBounds(bounds).getCenter(), 5, { animate: false });
    map.invalidateSize({ pan: false });
    map.fitBounds(bounds, { animate: false });
  }
  function showTab(name) {
    tabs.forEach(function (a) { a.classList.toggle('active', a.dataset.tab === name); });
    document.getElementById('tab-forecast').hidden = name !== 'forecast';
    document.getElementById('tab-about').hidden = name !== 'about';
    if (name === 'forecast' && map) {
      if (pendingBounds) fitMap(pendingBounds);
      else map.invalidateSize();
    }
  }
  tabs.forEach(function (a) {
    a.addEventListener('click', function (e) {
      e.preventDefault();
      history.replaceState(null, '', '#' + a.dataset.tab);
      showTab(a.dataset.tab);
    });
  });
  showTab(location.hash === '#about' ? 'about' : 'forecast');

  // ---- Map shell (tiles render while the data loads) ----
  map = L.map('map', { preferCanvas: true, zoomSnap: 0.25 });
  fitMap([[-20.27, 120.56], [-15.44, 124.86]]); // replaced by the data's bounds

  L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Topo_Map/MapServer/tile/{z}/{y}/{x}', {
    opacity: 0.75,
    attribution: 'Tiles &copy; Esri &mdash; Esri, DeLorme, NAVTEQ, TomTom, Intermap, iPC, USGS, FAO, ' +
      'NPS, NRCAN, GeoBase, Kadaster NL, Ordnance Survey, Esri Japan, METI, Esri China (Hong Kong), ' +
      'and the GIS User Community'
  }).addTo(map);

  // Points sit above the boundary polygons so polygons never swallow their clicks.
  map.createPane('points');
  map.getPane('points').style.zIndex = 450;
  var pointRenderer = L.canvas({ pane: 'points', padding: 0.5 });

  // Legend
  var legend = L.control({ position: 'bottomright' });
  legend.onAdd = function () {
    var div = L.DomUtil.create('div', 'legend');
    var html = '<strong>Probability of<br>toad arrival</strong>';
    for (var i = 0; i < COLOURS.length; i++) {
      html += '<div class="row-x"><i style="background:' + COLOURS[i] + '"></i>' +
        BREAKS[i].toFixed(1) + ' &ndash; ' + BREAKS[i + 1].toFixed(1) + '</div>';
    }
    div.innerHTML = html;
    return div;
  };
  legend.addTo(map);

  // Fullscreen button (stands in for bslib's card full_screen = TRUE)
  var card = document.querySelector('.forecast-map');
  if (card.requestFullscreen) {
    var fs = L.control({ position: 'topleft' });
    fs.onAdd = function () {
      var div = L.DomUtil.create('div', 'leaflet-bar');
      var a = L.DomUtil.create('a', 'fullscreen-btn', div);
      a.href = '#'; a.title = 'Full screen'; a.setAttribute('role', 'button');
      a.innerHTML = '&#x26F6;';
      L.DomEvent.on(a, 'click', function (e) {
        L.DomEvent.preventDefault(e);
        if (document.fullscreenElement) document.exitFullscreen();
        else card.requestFullscreen();
      });
      L.DomEvent.disableClickPropagation(div);
      return div;
    };
    fs.addTo(map);
    document.addEventListener('fullscreenchange', function () { map.invalidateSize(); });
  }

  // ---- Load data ----
  function getJSON(url) {
    return fetch(url).then(function (r) {
      if (!r.ok) throw new Error(url + ': HTTP ' + r.status);
      return r.json();
    });
  }

  var status = document.getElementById('map-status');

  getJSON('data/forecast.json').then(function (fc) {
    init(fc);
    status.hidden = true;
    // Boundaries are context only; load them after the points are on screen.
    return Promise.all([getJSON('data/tcz.geojson'), getJSON('data/pastoral.geojson')]);
  }).then(function (geo) {
    addBoundaries(geo[0], geo[1]);
  }).catch(function (err) {
    console.error(err);
    status.hidden = false;
    status.textContent = 'Sorry, the forecast data could not be loaded. Please try refreshing the page.';
  });

  function addBoundaries(tcz, pastoral) {
    var tczLayer = L.geoJSON(tcz, {
      style: { color: '#cc0000', weight: 2, fillOpacity: 0.05 }
    }).bindTooltip('Toad Containment Zone', { sticky: true });
    var pastoralLayer = L.geoJSON(pastoral, {
      style: { color: '#1a6faf', weight: 1, fillOpacity: 0.05 }
    });
    // Same draw order as the R app: TCZ first, pastoral on top.
    tczLayer.addTo(map);
    pastoralLayer.addTo(map);
    L.control.layers(null, {
      'TCZ Boundary': tczLayer,
      'Pastoral Boundaries': pastoralLayer
    }, { collapsed: false }).addTo(map);
  }

  function init(fc) {
    var years = fc.years, pts = fc.points, n = pts.lon.length;

    document.querySelectorAll('.n-sims').forEach(function (el) { el.textContent = fc.n_sims; });
    document.querySelectorAll('.updated').forEach(function (el) { el.textContent = fc.updated_label; });
    fitMap(fc.bounds);

    // One circle marker per waterpoint, restyled per year rather than rebuilt.
    var markers = new Array(n);
    var yearIdx = 0;
    for (var i = 0; i < n; i++) {
      var m = L.circleMarker([pts.lat[i], pts.lon[i]], {
        renderer: pointRenderer, stroke: false, radius: GREY.radius,
        fillColor: GREY.fillColor, fillOpacity: GREY.fillOpacity
      });
      m._idx = i;
      m.bindPopup(popupHtml);
      markers[i] = m;
    }
    var group = L.layerGroup(markers).addTo(map);

    function popupHtml(layer) {
      var i = layer._idx, yr = years[yearIdx], p = fc.prob[yearIdx][i];
      var html = '<b>Waterpoint</b><br>Probability colonised by <b>' + yr + '</b>: <b>' +
        Math.round(p * 100) + '%</b>';
      if (p > 0) {
        var arr = pts.arrival[i];
        html += '<br>' + (arr <= fc.start_year
          ? 'Colonised since: <b>' + arr + '</b>'
          : 'Median expected arrival year: <b>' + arr + '</b>');
      }
      return html;
    }

    function draw(idx) {
      yearIdx = idx;
      var probs = fc.prob[idx];
      for (var i = 0; i < n; i++) {
        var p = probs[i];
        markers[i].setStyle(p > 0
          ? { radius: 4, fillColor: colourFor(p), fillOpacity: 0.85 }
          : GREY);
      }
      // Coloured dots on top of grey ones.
      for (var j = 0; j < n; j++) if (probs[j] > 0) markers[j].bringToFront();
      // Keep an open popup in step with the year.
      group.eachLayer(function (m) { if (m.isPopupOpen()) m.setPopupContent(popupHtml(m)); });
      readout.textContent = years[idx];
      slider.setAttribute('aria-valuetext', years[idx]);
    }

    // ---- Year slider + play ----
    var slider = document.getElementById('year');
    var readout = document.getElementById('year-readout');
    slider.min = 0; slider.max = years.length - 1; slider.value = 0;
    document.getElementById('year-min').textContent = years[0];
    document.getElementById('year-max').textContent = years[years.length - 1];
    slider.addEventListener('input', function () { stop(); draw(+slider.value); });

    var playBtn = document.getElementById('play'), timer = null;
    function setPlayButton(playing) {
      playBtn.querySelector('.play-icon').innerHTML = playing ? '&#10074;&#10074;' : '&#9654;';
      playBtn.querySelector('.play-label').textContent = playing ? 'Pause' : 'Play';
      playBtn.setAttribute('aria-label', playing ? 'Pause animation' : 'Play animation');
    }
    function stop() {
      if (timer) { clearInterval(timer); timer = null; setPlayButton(false); }
    }
    function step() {
      var next = +slider.value + 1;
      if (next > +slider.max) { stop(); return; }
      slider.value = next;
      draw(next);
      if (next === +slider.max) stop();
    }
    playBtn.addEventListener('click', function () {
      if (timer) { stop(); return; }
      if (+slider.value === +slider.max) { slider.value = 0; draw(0); }
      setPlayButton(true);
      timer = setInterval(step, PLAY_MS);
    });

    draw(0);

    // Handle for debugging from the browser console.
    window.tczForecast = { map: map, markers: markers, setYear: function (i) { slider.value = i; draw(i); } };
  }
})();
