import React from 'react';
import './Badge.style.scss';

const Badge = (props) => {
  const { children, isActive } = props;
  return <div className={`badge ${isActive ? 'active' : ''}`}>{children}</div>;
};

export default Badge;
