import React, { useState } from 'react';
import './ButtonToggle.style.scss';

const ButtonToggle = (props: {
  items: any[];
  value?: any;
  labelKey: string;
  onChange?: any;
  label: string;
}) => {
  const {
    items = [],
    value,
    labelKey = '',
    onChange = () => {},
    label = '',
  } = props;
  const [currentItem, setCurrentItem] = useState(value || items[0]);

  const handleButtonClick = (item: any) => {
    setCurrentItem(item);
    onChange(item);
  };

  return (
    <div className="form-element">
      <label className="form__label">{label}</label>
      <div className=" button-toggle">
        {items.map((item, index: number) => (
          <div
            key={index}
            className={`button-toggle__item ${
              item === currentItem ? 'selected' : ''
            } `}
            onClick={() => handleButtonClick(item)}
          >
            {item[labelKey]}
          </div>
        ))}
      </div>
    </div>
  );
};

export default ButtonToggle;
