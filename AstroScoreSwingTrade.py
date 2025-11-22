import streamlit as st
import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.coordinates import get_body, solar_system_ephemeris, EarthLocation
import astropy.units as u
from datetime import datetime
import yfinance as yf
import plotly.graph_objects as go
import warnings
warnings.filterwarnings('ignore')

# === SAME NATAL DATA AS BEFORE (BSE, NSE, Reliance, HDFC Bank, etc.) ===
NATAL_DATA = { ... }  # Keep the same dictionary from previous version

# === ALL PREVIOUS FUNCTIONS (get_planet_positions, aspect_score, etc.) ===
# Copy-paste all functions from the previous app here (unchanged)

# ==================== NEW: EPHEMERIS WHEEL FUNCTION ====================
def draw_astro_wheel(natal_pos, transits, title="BSE + Transits Wheel"):
    fig = go.Figure()

    # Background circles
    angles = np.linspace(0, 2*np.pi, 100)
    fig.add_trace(go.Scatter(x=np.cos(angles), y=np.sin(angles), mode='lines', line=dict(color='gray', width=3), name='Zodiac'))
    fig.add_trace(go.Scatter(x=0.9*np.cos(angles), y=0.9*np.sin(angles), mode='lines', line=dict(color='lightgray'), name='Inner'))

    # Signs labels
    signs = ["Aries", "Taurus", "Gemini", "Cancer", "Leo", "Virgo", "Libra", "Scorpio", "Sagittarius", "Capricorn", "Aquarius", "Pisces"]
    for i, sign in enumerate(signs):
        angle = (i * 30 + 15) * np.pi / 180
        fig.add_annotation(x=1.15*np.cos(angle), y=1.15*np.sin(angle), text=sign[:3], showarrow=False, font=dict(size=10))

    # Planets symbols
    symbols = {
        'sun': '☉', 'moon': '☽', 'mercury': '☿', 'venus': '♀', 'mars': '♂',
        'jupiter': '♃', 'saturn': '♄', 'uranus': '♅', 'neptune': '♆', 'pluto': '♇'
    }

    # Plot NATAL planets (blue)
    for name, lon in natal_pos.items():
        if name in ['Asc', 'Sun', 'Moon', 'Mars', 'Jupiter']:
            rad = lon * np.pi / 180
            x, y = np.cos(rad), np.sin(rad)
            fig.add_trace(go.Scatter(x=[x], y=[y], mode='text+markers',
                                   marker=dict(size=18, color='blue', symbol='circle'),
                                   text=symbols.get(name.lower(), name[0]), textposition="middle center",
                                   name=f"Natal {name}", textfont=dict(size=16)))

    # Plot TRANSIT planets (red)
    for planet, lon in transits.items():
        rad = lon * np.pi / 180
        x, y = 0.85 * np.cos(rad), 0.85 * np.sin(rad)  # Slightly inside
        fig.add_trace(go.Scatter(x=[x], y=[y], mode='text+markers',
                               marker=dict(size=16, color='red', symbol='diamond'),
                               text=symbols.get(planet, planet[0].upper()),
                               name=f"{planet.capitalize()} Transit", textfont=dict(size=14)))

    # Draw aspects (only strong ones ≤8°)
    aspect_colors = {0: 'green', 90: 'red', 180: 'red', 120: 'green', 60: 'green'}
    for natal_name, n_lon in natal_pos.items():
        for p, t_lon in transits.items():
            diff = min(abs(n_lon - t_lon), 360 - abs(n_lon - t_lon))
            if diff <= 8:
                aspect_type = min([a for a in [0,60,90,120,180] if abs(diff - a) <= 8], default=None)
                if aspect_type:
                    rad1 = n_lon * np.pi/180
                    rad2 = t_lon * np.pi/180
                    fig.add_trace(go.Scatter(x=[np.cos(rad1), np.cos(rad2)*0.85],
                                           y=[np.sin(rad1), np.sin(rad2)*0.85],
                                           mode='lines',
                                           line=dict(color=aspect_colors.get(aspect_type, 'gray'), dash='dot' if aspect_type in [60,120] else 'solid'),
                                           name=f"{p}-{natal_name} {aspect_type}°"))

    fig.update_layout(
        title=title,
        width=700, height=700,
        xaxis=dict(range=[-1.3,1.3], showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(range=[-1.3,1.3], showgrid=False, zeroline=False, showticklabels=False),
        showlegend=True,
        paper_bgcolor='black',
        plot_bgcolor='black',
        font=dict(color='white')
    )
    return fig

# ==================== MAIN APP WITH WHEEL ====================
st.title("🔭 Jensen Astro-Cycles Pro App + Ephemeris Wheel")
st.markdown("**Real-time 20–30 day Nifty/BSE forecast + Live Astrology Wheel**")

# Sidebar same as before...
# ... (market selection code)

if st.button("🚀 Calculate + Show Wheel"):
    # All previous calculations
    loc = {'lat': NATAL_DATA[market]['lat'], 'lon': NATAL_DATA[market]['lon']}
    t_now = Time.now()
    transits = get_planet_positions(t_now, loc)
    natal_pos = compute_natal_positions(market)
    score, aspects = aspect_score(natal_pos, transits)
    triggers = check_triggers(transits)
    levels = get_price_levels("BSESN" if "BSE" in market else "^NSEI")
    action, notes = suggest_action(score, triggers, levels)

    # === DISPLAY RESULTS ===
    col1, col2 = st.columns([1, 1])

    with col1:
        st.subheader(f"🔮 20-30 Day Trend – {market}")
        st.success(action)
        st.info(notes)
        st.metric("Jensen Score", f"{score:+.1f}")
        if triggers:
            st.warning("⚡ Triggers: " + " | ".join(triggers))

    with col2:
        st.subheader("🌌 Live Ephemeris Wheel")
        wheel = draw_astro_wheel(natal_pos, transits, f"{market} – {datetime.now().strftime('%d %b %Y %H:%M')} IST")
        st.plotly_chart(wheel, use_container_width=True)

    # Tables
    st.subheader("Key Applying Aspects")
    for a in aspects[:15]:
        st.write("• " + a)