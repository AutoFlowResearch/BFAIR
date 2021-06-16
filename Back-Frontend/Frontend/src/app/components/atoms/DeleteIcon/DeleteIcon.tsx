import React from 'react';
import deleteIconImg from '../../../../assets/images/delete-icon.svg';

const DeleteIcon = (props) => {
  const { onClick = () => {} } = props;

  const handleClick = (event) => {
    event && event.stopPropagation && event.stopPropagation();
    onClick();
  };

  return (
    <img
      alt=''
      src={deleteIconImg}
      onClick={handleClick}
    />
  );
};
export default DeleteIcon;
